# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import sys
import h5py
import numpy as np
import contextlib
from datetime import datetime

from . import h5fs

class NexusException(Exception):
    """
    Base class for generic Nexus exceptions.
    """
    pass

class NexusFormatException(NexusException):
    """
    Raised when the hdf5 structure doesn't match the Nexus standards.
    """
    pass

if sys.version_info < (3,):
    text_dtype = h5py.special_dtype(vlen=unicode)
else:
    text_dtype = h5py.special_dtype(vlen=str)

def astextarray(arr):
    return np.asarray(arr, dtype=text_dtype)

def timestamp():
    return datetime.now().isoformat()
    
    
class Path(h5fs.Path):

    @property
    def nxclass(self):
        return self.stats().get('NX_class',None)
    
    def updated(self):
        self.nxroot().update_stats(file_update_time=timestamp())
    
    def _init_nxclass(self,path,nxclass,statsgen=None,**openargs):
        with self.h5open(**openargs):
            path = path.mkdir()
            with path.open() as node:
                nxclassattr = node.attrs.get('NX_class',None)
                
                if nxclassattr:
                    if nxclassattr!=nxclass:
                        raise NexusFormatException('NX_class=={} (should be {})'
                                                   .format(nxclassattr,nxclass))
                else:
                    node.attrs["NX_class"] = nxclass
                    if statsgen:
                        node.attrs.update(**statsgen())
                    path.updated()
            return path
    
    def _verify_class(self,*nxclasses):
        nxclass = self.nxclass
        if all(nxclass!=cls for cls in nxclasses):
            raise NexusFormatException('NX_class=={} while it should be one of these: {}'
                                       .format(nxclass,nxclasses))
    
    def _verify_not_class(self,*nxclasses):
        nxclass = self.nxclass
        if any(nxclass==cls for cls in nxclasses):
            raise NexusFormatException('NX_class=={} while it shouldn\'t be one of these: {}'
                                       .format(nxclass,nxclasses))
    
    def _verify_has_class(self):
        if not self.nxclass:
            raise NexusFormatException('NX_class attribute is missing')
    
    def _lookup(self,nxclass):
        with self.h5open(mode='r'):
            parent = self.parent
            if not parent:
                return self
            while parent.nxclass and parent.parent and \
                  parent.nxclass!=nxclass and parent!=parent.parent:
                parent = parent.parent
            return parent

    def _statsgen_nxroot(self):
        return {'file_name':self._h5file.path,
                'file_time':timestamp(),
                'HDF5_Version':h5py.version.hdf5_version,
                'h5py_version':h5py.version.version,
                'creator':'spectrocrunch'}

    def nxroot(self,**openargs):
        return self._init_nxclass(self.root,'NXroot',
                                  statsgen=self._statsgen_nxroot,
                                  **openargs)
    
    def nxentry(self,name=None,**openargs):
        path = self._lookup('NXentry')
        if path.nxclass=='NXentry':
            if not name:
                return path
            if path.name==name:
                return path
            path = path.parent
        elif name is None:
            raise NexusException('Could not find NXentry')
        path._verify_class('NXroot')
        return self._init_nxclass(path[name],'NXentry',**openargs)

    def nxsubentry(self,name,**openargs):
        self._verify_class('NXentry','NXsubentry')
        return self._init_nxclass(self[name],'NXsubentry',**openargs)
    
    def nxdata(self,name,**openargs):
        self._verify_not_class('NXroot')
        path = self._init_nxclass(self[name],'NXdata',**openargs)
        return self.factory(path,cls=_NXData)
        
    def nxinstrument(self,**openargs):
        path = self._lookup('NXinstrument')
        if path.nxclass=='NXinstrument':
            return path
        
        path = self._lookup('NXentry')
        if path.nxclass!='NXentry':
            raise NexusException('Could not find NXinstrument or NXentry')
        
        return self._init_nxclass(path['instrument'],'NXinstrument',**openargs)
     
    def factory(self,path,cls=None):
        if cls is None:
            return super(Path,self).factory(path)
        else:
            return cls(path,h5file=self._h5file)
        
        
class _NXData(Path):
    
    @contextlib.contextmanager
    def _verify(self):
        with self.h5open():
            self._verify_class('NXdata')
            yield
    
    @property
    def signal(self):
        with self._verify():
            with self.open() as node:
                name = node.attrs.get("signal",None)
                if name:
                    return self[name]
                else:
                    raise NexusFormatException('NXdata does not have a signal')
    
    def signals(self):
        yield self.signal
        for signal in self.auxiliary_signals:
            yield signal
            
    def auxiliary_signals(self):
        with self._verify():
            with self.open() as node:
                names = node.attrs.get("auxiliary_signals",[])
                for name in names:
                    yield self[name]
    
    def add_signal(self,name,title=None,**createparams):
        with self._verify():
            if title is None:
                title = name
            with self[name].open(**createparams) as node:
                 node.attrs["long_name"] = astextarray(title)
            
            with self.open() as node:
                signals = np.append(node.attrs.get("auxiliary_signals",[]),
                                    node.attrs.get("signal",[]))
                if signals.size:
                    node.attrs["auxiliary_signals"] = signals.astype(text_dtype)
                node.attrs["signal"] = astextarray(name)
 
    def default_signal(self,name):
        with self._verify():
            with self.open() as node:
                if name==node.attrs["signal"]:
                    return
                signals = np.append(node.attrs.get("auxiliary_signals",[]),
                                    node.attrs.get("signal",[]))
                if name not in signals:
                    raise ValueError('No signal with that name')
                                    
                aux = astextarray([signal for signal in signals if signal!=name])
                if aux.size:
                    node.attrs["auxiliary_signals"] = aux
                node.attrs["signal"] = astextarray(name)

    def get_signal(self,name):
        with self._verify():
            for signal in self.signals():
                if signal.name == name:
                    return signal
            raise fs.Missing(self[name])
    
    def remove_signal(self,name):
        with self._verify():
            self[name].remove()
            with self.open() as node:
                signals = node.attrs.pop("auxiliary_signals",None)
                if signals.size:
                    if self.signal.name==name:
                        aux = signals[:-1]
                        if aux.size:
                            node.attrs["auxiliary_signals"] = aux
                        node.attrs["signal"] = signals[-1]
                    else:
                        aux = astextarray([signal for signal in signals if signal!=name])
                        if aux.size:
                            node.attrs["auxiliary_signals"] = aux

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
import json

from . import h5fs
from ..common import units
from ..common import instance
from .. import __version__

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

def textarray(arr):
    return np.array(arr, dtype=text_dtype)

def timestamp():
    return textarray(datetime.now().isoformat())
    
    
class Path(h5fs.Path):

    def lsinfo(self):
        name,stats = super(Path,self).lsinfo()
        nxclass = stats.pop('NX_class',None)
        if nxclass:
            name = '{}:{}'.format(name,nxclass)
        return name,stats
        
    @property
    def nxclass(self):
        if self.exists:
            return self.stats().get('NX_class',None)
        else:
            return None
            
    def updated(self):
        with self.h5open():
            tm = timestamp()
            for path in self.iterup_is_nxclass('NXentry','NXsubentry'):
                path['end_time'].mkfile(mode='w',data=tm)
            self.nxroot().update_stats(file_update_time=tm)
    
    def _init_nxclass(self,path,nxclass,attrgen=None,filesgen=None,**openparams):
        with self.h5open(**openparams):
            path = path.mkdir()
            with path.open() as node:
                nxclassattr = node.attrs.get('NX_class',None)
                
                if nxclassattr:
                    if nxclassattr!=nxclass:
                        raise NexusFormatException('NX_class=={} (should be {})'
                                                   .format(nxclassattr,nxclass))
                else:
                    node.attrs["NX_class"] = nxclass
                    if attrgen:
                        node.attrs.update(**attrgen())
                    if filesgen:
                        for name,data in filesgen().items():
                            path[name].mkfile(data=data,**openparams)
                            
                    path.updated()
            return self.factory(path)
    
    def is_nxclass(self,*nxclasses):
        nxclass = self.nxclass
        if nxclasses:
            return any(nxclass==cls for cls in nxclasses)
        else:
            return bool(nxclass)
    
    def is_not_nxclass(self,*nxclasses):
        nxclass = self.nxclass
        if nxclasses:
            return all(nxclass!=cls for cls in nxclasses)
        else:
            return not bool(nxclass)
            
    def _raise_if_class(self,*nxclasses):
        if self.is_nxclass(*nxclasses):
            raise NexusFormatException('NX_class should be one of these: {}'
                                       .format(nxclasses))
    
    def _raise_ifnot_class(self,*nxclasses):
        if self.is_not_nxclass(*nxclasses):
            raise NexusFormatException('NX_class shouldn\'t be one of these: {}'
                                       .format(nxclasses))
    
    def iterup_is_nxclass(self,*nxclasses):
        with self.h5open(mode='r'):
            for path in self.iterup:
                if path.is_nxclass(*nxclasses):
                    yield path

    def iter_is_nxclass(self,*nxclasses):
        with self.h5open(mode='r'):
            for path in self:
                if path.is_nxclass(*nxclasses):
                    yield path

    def iterup_isnot_nxclass(self,*nxclasses):
        with self.h5open(mode='r'):
            for path in self.iterup:
                if path.is_not_nxclass(*nxclasses):
                    yield path

    def findfirstup_is_nxclass(self,*nxclasses):
        path = None
        for path in self.iterup_is_nxclass(*nxclasses):
            return path
        return self.nxroot()

    def nxroot(self,**openparams):
        return self._init_nxclass(self.root,'NXroot',
                                  attrgen=self._attrgen_nxroot,
                                  **openparams)
    
    def _attrgen_nxroot(self):
        return {'file_name':textarray(self._h5file.path),
                'file_time':timestamp(),
                'HDF5_Version':textarray(h5py.version.hdf5_version),
                'h5py_version':textarray(h5py.version.version),
                'creator':textarray('spectrocrunch')}
                
    def nxentry(self,name=None,**openparams):
        path = self.findfirstup_is_nxclass('NXentry')
        if path.nxclass=='NXentry':
            if not name:
                return path
            if path.name==name:
                return path
            path = path.parent
        elif name is None:
            raise NexusException('Could not find NXentry')
        path._raise_ifnot_class('NXroot')
        return self._init_nxclass(path[name],'NXentry',
                                  filesgen=self._files_nxentry,
                                  **openparams)

    def nxsubentry(self,name,**openparams):
        self._raise_ifnot_class('NXentry','NXsubentry')
        return self._init_nxclass(self[name],'NXsubentry',
                                  filesgen=self._files_nxentry,
                                  **openparams)
    
    def _files_nxentry(self):
        return {'start_time':timestamp()}

    def nxdata(self,name,**openparams):
        self._raise_if_class('NXroot')
        return self._init_nxclass(self[name],'NXdata',**openparams)
    
    def nxnote(self,name,**openparams):
        self._raise_if_class('NXroot')
        return self._init_nxclass(self[name],'NXnote',
                                  attrgen=self._attrgen_nxnote,
                                  **openparams)
    
    def _attrgen_nxnote(self):
        return {'date':timestamp()}
        
    def nxprocess(self,name,**openparams):
        self._raise_ifnot_class('NXentry')
        return self._init_nxclass(self[name],'NXprocess',
                                  filesgen=self._filesgen_nxprocess,
                                  **openparams)
    
    def _filesgen_nxprocess(self):
        i = 0
        for path in self.iter_is_nxclass('NXprocess'):
            ind = path['sequence_index']
            if ind.exists:
                with ind.open(mode='r') as node:
                    i = max(i,node.value)
        return {'program':textarray('spectrocrunch'),
                'version':textarray(__version__),
                'date':timestamp(),
                'sequence_index':i+1}

    def nxinstrument(self,**openparams):
        path = self.findfirstup_is_nxclass('NXinstrument')
        if path.nxclass=='NXinstrument':
            return path
        
        path = self.findfirstup_is_nxclass('NXentry')
        if path.nxclass!='NXentry':
            raise NexusException('Could not find NXinstrument or NXentry')
        
        return self._init_nxclass(path['instrument'],'NXinstrument',**openparams)
    
    def nxcollection(self,name,**openparams):
        self._raise_if_class('NXroot')
        return self._init_nxclass(self[name],'NXcollection',**openparams)

    def positioners(self,**openparams):
        path = self.nxinstrument(**openparams)
        return path.nxcollection('positioners',**openparams)
        
    def factory(self,path):
        cls = None
        if isinstance(path,h5fs.Path):
            nxclass = path.nxclass
            if nxclass=='NXdata':
                cls = _NXdata
            elif nxclass=='NXprocess':
                cls = _NXprocess
            elif nxclass=='NXnote':
                cls = _NXnote
            elif nxclass=='NXcollection':
                if path.name == 'positioners' and path.parent.nxclass == 'NXinstrument':
                    cls = _Positioners
        
        if cls is None:
            return super(Path,self).factory(path)
        else:
            return cls(path,h5file=self._h5file)
    
    def mark_default(self):
        path = self
        parent = path.parent
        while parent.nxclass:
            parentnxclass = parent.nxclass
            if parentnxclass=='NXdata':
                parent.default_signal(path.name)
            else:
                parent.update_stats(default=path.name)
            if parentnxclass=='NXEntry' or parentnxclass=='NXroot':
                break
            path = parent
            parent = path.parent
    
    def default(self):
        path = self.nxentry()
        default = entry.stats().get('default',None)
        while default:
            path = path['default']
            default = entry.stats().get('default',None)
        return path
        
        
class _NXprocess(Path):
    
    @property
    def config(self):
        return self.nxnote('configuration')
    
    @property
    def results(self):
        return self.nxcollection('results')

    @property
    def plotselect(self):
        return self.nxdata('plotselect')
        
        
class _NXnote(Path):
    
    @property
    def data(self):
        with self['data'].open(mode='r') as node:
            data = node.value
            mimetype = node.attrs.get('type','text/plain')
            if mimetype=='application/json':
                data = json.loads(data)
            return data
            
    def _write(self,data,mimetype):
        with self['data'].open(data=data) as node:
            node.attrs['type'] = textarray(mimetype)
            
    def write_text(self,text):
        self._write(textarray(text),'text/plain')
    
    def write_json(self,dic):
        self._write(textarray(json.dumps(dic)),'application/json')
        
        
class _NXdata(Path):
    
    @contextlib.contextmanager
    def _verify(self):
        with self.h5open():
            self._raise_ifnot_class('NXdata')
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
    
    @property
    def signals(self):
        yield self.signal
        for signal in self.auxiliary_signals:
            yield signal
    
    @property
    def auxiliary_signals(self):
        with self._verify():
            with self.open() as node:
                names = node.attrs.get("auxiliary_signals",[])
                for name in names:
                    yield self[name]
    
    def add_signal(self,name=None,path=None,title=None,interpretation=None,**createparams):
        """
        args:
            name(str)
            signal(Path)
            title(str)
            createparams(dict)
        """
        with self._verify():
            # Create the signal dataset
            if path:
                if name is None:
                    name = path.name
                self[name].link(path)
            elif name:
                if title is None:
                    title = name
                with self[name].open(**createparams) as node:
                    node.attrs["long_name"] = textarray(title)
                    if interpretation:
                        node.attrs["interpretation"] = textarray(interpretation)
            else:
                raise ValueError('Provide either "name" or "signal"')
            
            # Add signal name to NXdata attributes
            with self.open() as node:
                signals = np.append(node.attrs.get("auxiliary_signals",[]),
                                    node.attrs.get("signal",[]))
                if signals.size:
                    node.attrs["auxiliary_signals"] = signals.astype(text_dtype)
                node.attrs["signal"] = textarray(name)
            
            self.updated()
            return self[name]
    
    def default_signal(self,name):
        with self._verify():
            with self.open() as node:
                if name==node.attrs["signal"]:
                    return
                signals = np.append(node.attrs.get("auxiliary_signals",[]),
                                    node.attrs.get("signal",[]))
                if name not in signals:
                    raise ValueError('No signal with that name')
                                    
                aux = textarray([signal for signal in signals if signal!=name])
                if aux.size:
                    node.attrs["auxiliary_signals"] = aux
                node.attrs["signal"] = textarray(name)
    
            self.updated()
                
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
                        aux = textarray([signal for signal in signals if signal!=name])
                        if aux.size:
                            node.attrs["auxiliary_signals"] = aux
            self.updated()

    def set_axes(self,*axes):
        """
        Args:
            *axes(3-tuple|str): (name(str),value(array|Path|None),attr(dict|None)) or name(str)
        """
        with self._verify():
            with self.open() as node:
                # Look for the signal
                signal = node.attrs.get('signal',None)
                if signal:
                    signal = node.get(signal,default=None)
                if not signal:
                    raise NexusFormatException('NXdata should have a signal before adding axes')
            
                # Check dimensions
                axes,positioners = self._axes_validate(axes,signal)
                
                # Create axes if needed
                names = []
                for name,value,attr in axes:
                    axis = positioners.add_axis(name,value,**attr)
                    if self[name].linkdest()!=axis:
                        self[name].link(axis)
                    names.append(name)
                    
                node.attrs["axes"] = textarray(names)
                
                self.updated()
                return axes
                
    def _axes_validate(self,axesin,signal):
        positioners = self.positioners()
        
        # Axes as list of tuples
        axes = []
        for axis in axesin:
            if instance.isstring(axis):
                name,value,attr = axis,None,{}
            else:
                name,value,attr = axis
                if attr is None:
                    attr = {}
            axes.append((name,value,attr))

        # Check compatibility with data dimensions
        shaped = signal.shape
        shapea = tuple(self._axis_size(name,value,positioners) for name,value,_ in axes)
        if shaped!=shapea:
            raise NexusFormatException('Wrong axes dimensions: {} instead of {}'.shape(shaped,shapea))
        
        return axes,positioners

    def _axis_size(self,name,value,positioners):
        if isinstance(value,h5fs.Path):
            with value.open(mode='r') as anode:
                return len(anode)
        elif value:
            return len(value)
        else:
            return len(positioners.get_axis(name))
            
    @property
    def axes(self):
        with self._verify():
            with self.open() as node:
                ret = []
                for name in node.attrs["axes"]:
                    with self[name].open(mode='r') as axis:
                        ret.append(self._axis_astuple(name,axis))
                return ret

    def _axis_astuple(self,name,axis):
        attrs = {'title':axis.attrs.get('long_name',None),
                 'units':axis.attrs.get('units',None)}
        values = axis.value
        values = units.Quantity(values,attrs['units'])
        return name,values,attrs
        

class _Positioners(Path):
    
    def add_axis(self,name,value,title=None,units=None):
        axis = self[name]
        if axis.exists:
            if value is not None:
                if not self._axis_equal(axis,value):
                    raise ValueError('Axis {} already exists'.format(name))
        else:
            if value is None:
                raise ValueError('Axis {} does not exist yet and needs data'.format(name))
            if isinstance(value,h5fs.Path):
                axis.link(value)
            else:
                if title is None:
                    title = name
                
                if instance.isquantity(value):
                    units = str(value.units)
                    value = value.magnitude

                with axis.open(data=value) as node:
                    node.attrs["long_name"] = textarray(title)
                    if units:
                        node.attrs["units"] = textarray(units)
            self.updated()
        return axis
    
    def get_axis(self,name):
        axis = self[name]
        if axis.exists:
            with axis.open(mode='r') as node:
                values = node.value
                u = node.attrs.get('units','dimensionless')
                values = units.Quantity(values,u)
            return values
            
    def _axis_equal(self,axis1,axis2):
        if isinstance(axis1,h5fs.Path):
            with axis1.open(mode='r') as node1:
                if isinstance(axis2,h5fs.Path):
                    with axis2.open(mode='r') as node2:
                        return np.array_equal(node1.value,node2.value)
                else:
                    return np.array_equal(node1.value,axis2)
        else:
            if isinstance(axis2,h5fs.Path):
                with axis2.open(mode='r') as node2:
                    return np.array_equal(axis1,node2.value)
            else:
                return np.array_equal(axis1,axis2)

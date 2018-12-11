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
import re
import logging

from . import fs
from . import h5fs
from ..utils import units
from ..utils import instance
from ..utils import hashing
from .. import __version__
PROGRAM_NAME = 'spectrocrunch'

logger = logging.getLogger(__name__)

class NexusException(h5fs.FileSystemException):
    """
    Base class for generic Nexus exceptions.
    """
    pass

class NexusFormatException(NexusException):
    """
    Raised when the hdf5 structure doesn't match the Nexus standards.
    """
    pass

class NexusProcessWrongHash(NexusException):
    """
    Raised when the NXprocess exists with a different hash
    """
    pass

if sys.version_info < (3,):
    text_dtype = h5py.special_dtype(vlen=unicode)
else:
    text_dtype = h5py.special_dtype(vlen=str)

def textarray(arr):
    return np.array(arr, dtype=text_dtype)

def dataprepare(data):
    if instance.isstring(data):
        data = textarray(data)
    return data
    
def timestamp():
    return textarray(datetime.now().isoformat())

def calc_checksum(dependencies,confighash):
    if dependencies:
        hashes = [prev.checksum for prev in dependencies]
    else:
        hashes = []
    hashes.append(confighash)
    return hashing.mergejhash(*hashes)

def nxprocess_extract(name):
    start,end,num = 0,0,'0'
    for m in re.finditer("[0-9]+",name):
        start,end,num = m.start(), m.end(), m.group()
    return start,end,int(num),start==end

def nxprocess_fmtname(name):
    start,end,num,nonumber = nxprocess_extract(name)
    if nonumber:
        fmt = name+'.{:d}'
    else:
        fmt = name[:start]+'{{:0{}d}}'.format(end-start)+name[end:]
    return fmt,num,nonumber

def nxprocess_regex(name):
    name = re.escape(name)
    start,end,num,nonumber = nxprocess_extract(name)
    if nonumber:
        pattern = '^'+name+'(\.[0-9]+)?$'
    else:
        pattern = '^'+name[:start]+'[0-9]+'+name[end:]+'$'
    return pattern

def incrementname(name):
    fmt,num,_ = nxprocess_fmtname(name)
    return fmt.format(num+1)

class Path(h5fs.Path):
       
    def mkfile(self,data=None,units=None,attributes=None,**params):
        if instance.isquantity(data):
            units = '{:~}'.format(data.units)
            data = data.magnitude
        if data is not None:
            params['data'] = dataprepare(data)
        with self.open(**params) as dset:
            if units is not None:
                dset.attrs['units'] = dataprepare(units)
            if attributes is not None:
                for k,v in attributes.items():
                    dset.attrs[k] = dataprepare(v)
        return self
    
    def update_stats(self,**stats):
        stats = {k:dataprepare(v) for k,v in stats.items()}
        super(Path,self).update_stats(**stats)
    
    def factory(self,path):
        if not isinstance(path,h5fs.Path):
            path = super(Path,self).factory(path)
            
        cls = None
        nxclass = path.nxclass
        if nxclass=='NXdata':
            cls = _NXdata
        elif nxclass=='NXprocess':
            cls = _NXprocess
        elif nxclass=='NXnote':
            cls = _NXnote
        elif nxclass=='NXmonochromator':
            cls = _NXmonochromator
        elif nxclass=='NXcollection':
            if path.name == 'positioners' and path.parent.nxclass == 'NXinstrument':
                cls = _Positioners
            else:
                cls = _NXcollection

        if cls is None:
            return path
        else:
            return cls(path,h5file=self._h5file)
            
    @property
    def nxclass(self):
        if self.exists:
            return self.get_stat('NX_class',default=None)
        else:
            return None
            
    def updated(self):
        with self.h5open():
            tm = timestamp()
            for path in self.iterup_is_nxclass('NXentry','NXsubentry'):
                path['end_time'].write(mode='w',data=tm)
            for path in self.iterup_is_nxclass('NXprocess'):
                path['date'].write(mode='w',data=tm)
            self.nxroot().update_stats(file_update_time=tm)
    
    def _init_nxclass(self,path,nxclass,nxattributes=None,nxfiles=None,**openparams):
        with self.h5open(**openparams):
            path = path.mkdir()
            with path.open() as node:
                nxclassattr = node.attrs.get('NX_class',None)
                if nxclassattr:
                    if nxclassattr!=nxclass:
                        raise NexusFormatException('{} already exists with the wrong class ({} instead of {})'
                                                   .format(path.name,nxclassattr,nxclass))
                else:
                    node.attrs['NX_class'] = nxclass
                    if nxattributes:
                        if instance.iscallable(nxattributes):
                            nxattributes = nxattributes()
                        nxattributes.pop('NX_class',None)
                        node.attrs.update(**nxattributes)
                    if nxfiles:
                        if instance.iscallable(nxfiles):
                            nxfiles = nxfiles()
                        for name,data in nxfiles.items():
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
            raise NexusFormatException('NX_class shouldn\'t be one of these: {}'
                                       .format(nxclasses))
            
    def _raise_ifnot_class(self,*nxclasses):
        if self.is_not_nxclass(*nxclasses):
            raise NexusFormatException('NX_class should be one of these: {}'
                                       .format(nxclasses))
    
    def iterup_is_nxclass(self,*nxclasses):
        with self.h5open():
            for path in self.iterup:
                if path.is_nxclass(*nxclasses):
                    yield self.factory(path)

    def iterup_isnot_nxclass(self,*nxclasses):
        with self.h5open():
            for path in self.iterup:
                if path.is_not_nxclass(*nxclasses):
                    yield self.factory(path)

    def iter_is_nxclass(self,*nxclasses):
        with self.h5open():
            for path in self:
                if path.is_nxclass(*nxclasses):
                    yield self.factory(path)

    def iter_nxentry(self):
        for entry in self.root.iter_is_nxclass('NXentry'):
            yield entry

    def iter_nxprocess(self,allentries=True):
        if allentries:
            for entry in self.iter_nxentry():
                for process in entry.iter_is_nxclass('NXprocess'):
                    yield process
        else:
            entry = self.findfirstup_is_nxclass('NXentry')
            if entry.is_nxclass('NXentry'):
                for process in entry.iter_is_nxclass('NXprocess'):
                    yield process
            
    def findfirstup_is_nxclass(self,*nxclasses):
        path = None
        for path in self.iterup_is_nxclass(*nxclasses):
            return path
        return self.root

    def findfirstup_is_name(self,*names):
        path = None
        for path in self.iterup:
            if path.name in names:
                return path
        return self.root
        
    def nxroot(self,**openparams):
        return self._init_nxclass(self.root,'NXroot',
                                  nxattributes=self._nxattributes_nxroot,
                                  **openparams)
    
    def _nxattributes_nxroot(self):
        return {'file_name':textarray(self._h5file.path),
                'file_time':timestamp(),
                'HDF5_Version':textarray(h5py.version.hdf5_version),
                'h5py_version':textarray(h5py.version.version),
                'creator':textarray(PROGRAM_NAME)}
                
    def nxentry(self,name=None,**openparams):
        """Get NXentry (first looks in parents) and create if it does not exist
        """
        path = self.findfirstup_is_nxclass('NXentry')
        if path.is_nxclass('NXentry'):
            if path.name==name or not name:
                return path
        elif name is None:
            raise NexusException('Could not find NXentry')
        return self._init_nxclass(self.root[name],'NXentry',
                                  nxfiles=self._nxfiles_nxentry,
                                  **openparams)

    def get_nxentry(self,name=None,**openparams):
        path = self.findfirstup_is_nxclass('NXentry')
        if path.is_nxclass('NXentry'):
            if not name:
                return path
            if path.name==name:
                return path
            path = path.parent
        elif name is None:
            raise NexusException('Could not find NXentry')
        if path==self.root:
            path = path.nxroot()
        path._raise_ifnot_class('NXroot')

    def last_nxentry(self):
        entry = None
        for entry in self.root.iter_is_nxclass('NXentry'):
            pass
        if entry:
            return entry
        else:
            return self.root
            
    def new_nxentry(self,entry=None):
        name = self.next_nxentry_name(entry=entry)
        return self.nxentry(name=name)

    def next_nxentry_name(self,entry=None):
        # If argument not an NXentry: find it in the parents
        if entry is None:
            entry = self.findfirstup_is_nxclass('NXentry')
        elif entry.is_not_nxclass('NXentry'):
            entry = entry.findfirstup_is_nxclass('NXentry')
            if entry.is_not_nxclass('NXentry'):
                entry = self.findfirstup_is_nxclass('NXentry')

        # Default NXentry
        if entry.is_not_nxclass('NXentry'):
            entry = self.last_nxentry()

        # Next entry name
        if entry.is_nxclass('NXentry'):
            name = entry.name
            parent = entry.parent
            while parent[name].exists:
                name = incrementname(name)
        else:
            name = 'entry0001'
        
        return name
            
    def nxsubentry(self,name,**openparams):
        self._raise_ifnot_class('NXentry','NXsubentry')
        return self._init_nxclass(self[name],'NXsubentry',
                                  nxfiles=self._nxfiles_nxentry,
                                  **openparams)
    
    def _nxfiles_nxentry(self):
        return {'start_time':timestamp()}
        
    def nxdata(self,name,**openparams):
        self._raise_if_class('NXroot')
        return self._init_nxclass(self[name],'NXdata',**openparams)
    
    def nxnote(self,name,**openparams):
        self._raise_if_class('NXroot')
        return self._init_nxclass(self[name],'NXnote',
                                  nxfiles=self._nxfiles_nxnote,
                                  **openparams)
    
    def _nxfiles_nxnote(self):
        return {'date':timestamp()}
        
    def nxprocess(self,name,parameters=None,dependencies=None,allentries=True,**openparams):
        """Creates NXentry and NXprocess when needed.
        
        Args:
            name(str):
            parameters(dict):
            dependencies(list): a list of object with the "checksum" attribute
            openparams(Optional(dict)): overwrite device open parameters (if allowed)
        
        Returns:
            Path: exists and name may be incremented
        """
        process = self.get_nxprocess(name,parameters=parameters,dependencies=dependencies,allentries=allentries)
        if process.exists:
            return process
        else:
            # Create NXentry when needed
            entry = process.nxentry(name=process.parent.name)
            # Create NXprocess
            process = self._init_nxclass(process,'NXprocess',
                                         nxfiles=self._nxfiles_nxprocess,
                                         **openparams)
            process.set_config(parameters,dependencies=dependencies)
            return process
    
    def get_nxprocess(self,name=None,parameters=None,dependencies=None,allentries=True):
        """Get the path of an NXprocess, identified by its parameters and dependencies (and optionally the name)
        
        Args:
            name(str):
            parameters(dict): 
            dependencies(list): a list of object with the "checksum" attribute
        
        Returns:
            Path: may not exists yet and name may be incremented
        """
        if name:
            match = re.compile(nxprocess_regex(name)).match
        else:
            match = lambda x:True
        checksum = calc_checksum(dependencies,hashing.calcjhash(parameters))
        for process in self.iter_nxprocess(allentries=allentries):    
            if process.verify_checksum(checksum) and match(process.name):
                return process
        
        # Determine NXentry (may not exists yet)
        if self.exists:
            entry = self.findfirstup_is_nxclass('NXentry')
            if entry.is_not_nxclass('NXentry'):
                entry = entry.next_nxentry_name()
        else:
            if self.root==self.parent:
                entry = self
            else:
                entry = self.next_nxentry_name()

        # Make sure we have a unique process name
        if not name:
            name = incrementname('process')
        while entry[name].exists:
            name = incrementname(name)

        return entry[name]
    
    def _nxfiles_nxprocess(self):
        return {'program':textarray(PROGRAM_NAME),
                'version':textarray(__version__),
                'date':timestamp()}

    def nxinstrument(self,**openparams):
        path = self.findfirstup_is_nxclass('NXinstrument')
        if path.nxclass=='NXinstrument':
            return path
        path = self.findfirstup_is_nxclass('NXentry')
        if path.nxclass!='NXentry':
            raise NexusException('Could not find NXinstrument or NXentry')
        return self._init_nxclass(path['instrument'],'NXinstrument',**openparams)
    
    def measurement(self,**openparams):
        path = self.findfirstup_is_name('measurement')
        if path.name=='measurement':
            return path
        path = self.findfirstup_is_nxclass('NXentry')
        if path.nxclass!='NXentry':
            raise NexusException('Could not find measurement or NXentry')
        return self._init_nxclass(path['measurement'],'NXcollection',**openparams)
        
    def nxcollection(self,name,nxattributes=None,nxfiles=None,**openparams):
        self._raise_if_class('NXroot')
        return self._init_nxclass(self[name],'NXcollection',
               nxattributes=nxattributes,nxfiles=nxfiles,**openparams)

    def nxdetector(self,name,nxattributes=None,nxfiles=None,**openparams):
        instrument = self.nxinstrument(**openparams)
        return self._init_nxclass(instrument[name],'NXdetector',
               nxattributes=nxattributes,nxfiles=nxfiles,**openparams)

    def nxmonochromator(self,name='monochromator',nxattributes=None,nxfiles=None,**openparams):
        instrument = self.nxinstrument(**openparams)
        return self._init_nxclass(instrument[name],'NXmonochromator',
               nxattributes=nxattributes,nxfiles=nxfiles,**openparams)

    def application(self,name,definition=None,nxattributes=None,nxfiles=None,**openparams):
        entry = self.nxentry(**openparams)
        if definition is None:
            definition = name
        if nxfiles is None:
            nxfiles = {}
        nxfiles['definition'] = nxfiles.get('definition',definition)
        return self._init_nxclass(self[name],'NXsubentry',
               nxattributes=nxattributes,nxfiles=nxfiles,**openparams)

    def positioners(self,**openparams):
        path = self.findfirstup_is_nxclass('NXprocess')
        if path.nxclass=='NXprocess':
            path = path.results
        else:
            path = self.nxinstrument(**openparams)
        return path.nxcollection('positioners',**openparams)

    def mark_default(self):
        path = self
        parent = path.parent
        while parent.nxclass:
            parentnxclass = parent.nxclass
            if parentnxclass=='NXdata':
                parent.default_signal(path.name)
            elif parentnxclass in ['NXentry','NXsubentry','NXroot']:
                parent.update_stats(default=path.name)
            if parentnxclass=='NXroot':
                break
            path = parent
            parent = path.parent
            
    @property
    def default(self):
        path = self.root
        default = path.get_stat('default',default=None)
        while default:
            path = path[default]
            default = path.get_stat('default',default=None)
        return path
        
class _NXPath(Path):
    
    @contextlib.contextmanager
    def _verify(self):
        with self.h5open():
            self._verify_func()
            yield

    def _verify_func(self):
        return self._raise_ifnot_class(self.NX_CLASS)


class _NXprocess(_NXPath):
    
    NX_CLASS = 'NXprocess'
            
    @property
    def config(self):
        with self._verify():
            return self.nxnote('configuration')
    
    @property
    def configpath(self):
        return self['configuration']
        
    @property
    def results(self):
        with self._verify():
            return self.nxcollection('results')

    @property
    def resultspath(self):
        return self['results']
        
    @property
    def plotselect(self):
        with self._verify():
            return self.nxdata('plotselect')
    
    def set_config(self,parameters,dependencies=None):
        with self._verify():
            # Save parameters
            self.config.write_dict(parameters)

            # Links to dependencies
            if dependencies:
                #for dependency in dependencies:
                #    dependency._raise_ifnot_class(self.NX_CLASS)
                #    if self.parent!=dependency.parent:
                #        raise ValueError('{} and {} should be in the same entry'.format(self,dependencies))
                for dependency in dependencies:
                    dependency = self.dependencies[dependency.name].link(dependency)

            # Other info
            self['sequence_index'].write(data=self.sequence_index)
            self.configpath.update_stats(checksum=self.checksum)

            self.updated()
            
    def verify_checksum(self,checksum):
        return self.checksum==checksum
    
    @property
    def checksum(self):
        with self._verify():
            checksum = self.configpath.get_stat('checksum')
            if checksum is None:
                checksum = calc_checksum(self.dependencies,self.confighash)
            return checksum
                
    @property
    def confighash(self):
        if self.configpath.exists:
            return hashing.calcjhash(self.config.read(parse=True))
        else:
            return None
  
    @property
    def dependencies(self):
        return self.results.nxcollection('dependencies')
            
    @property
    def sequence_index(self):
        """Indexing starts from 1
        """
        path = self['sequence_index']
        if path.exists:
            return path.read()
        else:
            ind = 0
            with self._verify():
                for dependency in self.dependencies:
                    ind = max(ind,dependency.sequence_index)
            return ind+1
            
    
class _NXnote(_NXPath):
    
    NX_CLASS = 'NXnote'

    def read(self,parse=True):
        with self._verify():
            with self['data'].open(mode='r') as node:
                data = node[()]
                if parse:
                    mimetype = node.attrs.get('type','text/plain')
                    if mimetype=='application/json':
                        data = json.loads(data)
                return data
            
    def _write(self,data,mimetype):
        with self._verify():
            with self['data'].open(data=data) as node:
                node.attrs['type'] = textarray(mimetype)
            
    def write_text(self,text):
        with self._verify():
            self._write(textarray(text),'text/plain')
    
    def write_dict(self,dic):
        with self._verify():
            self._write(textarray(json.dumps(dic)),'application/json')

        
class _NXdata(_NXPath):
    
    NX_CLASS = 'NXdata'
    
    @property
    def signal(self):
        with self._verify():
            with self.open() as node:
                name = node.attrs.get("signal",None)
                if name:
                    return self[name]
                else:
                    return None
    
    @property
    def signal_shape(self):
        signal = self.signal
        if signal is None:
            return ()
        else:
            with signal.open(mode='r') as dset:
                return dset.shape
    
    @property
    def signals(self):
        with self._verify():
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
            # Check dimensions
            ashape = self.axes_shape
            if not ashape:
                ashape = self.signal_shape
            if ashape:
                if 'shape' in createparams:
                    dshape = tuple(createparams['shape'])
                elif 'data' in createparams:
                    dshape = createparams['data'].shape
                else:
                    if path is None:
                        spath = self[name]
                    else:
                        spath = path
                    with spath.open(mode='r') as dset:
                        dshape = dset.shape
                if dshape!=ashape:
                    raise NexusFormatException('Data dimensions {} do not match axes dimensions {}'
                                               .format(dshape,ashape))
            
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
                    if not interpretation:
                        interpretation = self._interpretation(node)
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
    
    def _interpretation(self,node):
        ndim = node.ndim
        if ndim==0:
            return 'scalar'
        elif ndim==1:
            return 'spectrum'
        else:
            return 'image'
    
    def default_signal(self,name):
        with self._verify():
            with self.open() as node:
                if name==node.attrs["signal"]:
                    return True
                signals = np.append(node.attrs.get("signal",[]),
                                    node.attrs.get("auxiliary_signals",[]))
                if name not in signals:
                    return False
                                    
                aux = textarray([signal for signal in signals if signal!=name])
                if aux.size:
                    node.attrs["auxiliary_signals"] = aux
                node.attrs["signal"] = textarray(name)
    
            self.updated()
        return True
        
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

    def sort_signals(self,other=None):
        with self._verify():
            with self.open() as node:
                signals = sorted(sig.name for sig in self.signals)
                
                if other:
                    other_signals = list(sig.name for sig in other.signals)
                    if sorted(other_signals)!=signals:
                        return
                    signals = other_signals

                aux = textarray(signals[1:])
                if aux.size:
                    node.attrs["auxiliary_signals"] = aux
                node.attrs["signal"] = textarray(signals[0])
    
            self.updated()

    def set_axes(self,*axes):
        """
        Args:
            *axes(3-tuple|str): (name(str),value(array|Path|None),attr(dict|None)) or name(str)
        """
        with self._verify():
            with self.open() as node:
                axes,positioners = self._axes_parse(axes)

                # Check dimensions
                dshape = self.signal_shape
                if dshape:
                    ashape = tuple(self._axis_size(name,value,positioners) for name,value,_ in axes)
                    if dshape!=ashape:
                        raise NexusFormatException('Axes dimensions {} do not match signal dimensions {}'
                                                   .format(ashape,dshape))

                # Create axes if needed
                names = []
                for name,value,attr in axes:
                    axis = positioners.add_axis(name,value,**attr)
                    if not self[name].exists:
                        self[name].link(axis,soft=True)
                    names.append(name)
                node.attrs["axes"] = textarray(names)
                self.updated()
                return axes

    def _axes_parse(self,axesin):
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
        
        return axes,positioners

    def _axis_size(self,name,value,positioners):
        if isinstance(value,h5fs.Path):
            with value.open(mode='r') as anode:
                return len(anode)
        elif value:
            return len(value)
        else:
            return positioners.axis_size(name)
            
    @property
    def axes(self):
        with self._verify():
            with self.open() as node:
                ret = []
                for name in node.attrs.get('axes',[]):
                    with self[name].open(mode='r') as axis:
                        ret.append(self._axis_astuple(name,axis))
                return ret

    @property
    def axes_shape(self):
        with self._verify():
            with self.open() as node:
                ret = []
                for name in node.attrs.get('axes',[]):
                    with self[name].open(mode='r') as axis:
                        ret.append(axis.size)
                return tuple(ret)

    def _axis_astuple(self,name,axis):
        attrs = {'title':axis.attrs.get('long_name',None),
                 'units':axis.attrs.get('units',None)}
        values = axis[()]
        values = units.Quantity(values,attrs['units'])
        return name,values,attrs
    
    def mark_default(self):
        super(_NXdata,self).mark_default()
        entry = self.nxentry()
        if self.parent!=entry:
            plotselect = entry['plotselect']
            if plotselect.islink:
                plotselect.remove()
            plotselect = plotselect.link(self)
            plotselect.mark_default()
            

class _NXcollection(_NXPath):
    
    NX_CLASS = 'NXcollection'
        
    def add_axis(self,name,value,title=None,units=None):
        with self._verify():
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
                    if instance.isquantity(value):
                        units = value.units
                        value = value.magnitude
                    if units:
                        try:
                            units = '{:~}'.format(units)
                        except ValueError:
                            pass
                    if title is None:
                        if units:
                            title = '{} ({})'.format(name,units)
                        else:
                            title = name
                    with axis.open(data=value) as node:
                        node.attrs["long_name"] = textarray(title)
                        if units:
                            node.attrs["units"] = textarray(units)
                self.updated()
            return axis
    
    def get_axis(self,name):
        with self._verify():
            axis = self[name]
            if axis.exists:
                with axis.open(mode='r') as node:
                    values = node[()]
                    u = node.attrs.get('units','dimensionless')
                    values = units.Quantity(values,u)
                return values
    
    def axis_size(self,name):
        with self._verify():
            axis = self[name]
            if axis.exists:
                with axis.open(mode='r') as node:
                    return node.size
            return None
            
    def _axis_equal(self,axis1,axis2):
        if isinstance(axis1,h5fs.Path):
            with axis1.open(mode='r') as node1:
                if isinstance(axis2,h5fs.Path):
                    with axis2.open(mode='r') as node2:
                        return np.array_equal(node1,node2)
                else:
                    return np.array_equal(node1,axis2)
        else:
            if isinstance(axis2,h5fs.Path):
                with axis2.open(mode='r') as node2:
                    return np.array_equal(axis1,node2)
            else:
                return np.array_equal(axis1,axis2)


class _Positioners(_NXcollection):
    
    NX_CLASS = 'NXcollection'
    
    def _verify_func(self):
        self._raise_ifnot_class(self.NX_CLASS)
        if self.name!='positioners':
            NexusFormatException('Name should be "positioners" not "{}"'
                                       .format(self.name))


class _NXmonochromator(_NXPath):
    
    NX_CLASS = 'NXmonochromator'
    
    @property
    def energy(self):
        with self._verify():
            with self['energy'].open(mode='r') as node:
                value = node[()]
                u = node.attrs.get('units','keV')
                return units.Quantity(value,u)
            
    @energy.setter
    def energy(self,value):
        with self._verify():
            self['energy'].mkfile(data=units.Quantity(value,'keV'))


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

import pandas as pd
from abc import ABCMeta, abstractmethod, abstractproperty
from future.utils import with_metaclass
import contextlib
import numpy as np

try:
    from contextlib import ExitStack
except ImportError:
    from contextlib2 import ExitStack
    
from . import axis
from ..utils import indexing
from ..io import nxfs
from ..math.interpolate import interpolate_regular

class RegularGrid(object):
    
    DEFAULT_STACKDIM = 0
    
    def __init__(self,axes,data,stackdim=None):
        self.axes = axes
        self.data = data
        if stackdim is None:
            self.stackdim = self.DEFAULT_STACKDIM
        else:
            self.stackdim = stackdim

    @property
    def shape(self):
        return tuple([ax.size for ax in self.axes])
    
    def __len__(self):
        return self.axes[0].size
        
    @property
    def ndim(self):
        return len(self.axes)

    @property
    def size(self):
        return np.prod(self.shape)
        
    def __repr__(self):
        return str(self.axes)
    
    def __str__(self):
        return str(self.axes)

    def sliceinfo(self,slicedim):
        """
        Args:
            slicedim(int)
        Returns:
            shape(tuple): shape after slicing
            indexgen(callable): slice index generator
        """
        if slicedim<0:
            slicedim += self.ndim
        maxdim = self.ndim-1

        if slicedim<0 or slicedim>maxdim:
            raise ValueError('Slice dimension should be between 0 and {}'.format(slicedim))

        if slicedim==0:
            indexgen = lambda i:(i,Ellipsis)
            shape = self.shape[1:]
            
        elif slicedim==maxdim:
            indexgen = lambda i:(Ellipsis,i)
            shape = self.shape[:-1]
        else:
            a = (slice(None),)*slicedim
            b = (slice(None),)*(maxdim-slicedim)
            indexgen = lambda i:a+(i,)+b
            shape = self.shape[:slicedim]+self.shape[slicedim+1:]
        
        return shape,indexgen
    
    @contextlib.contextmanager
    def open(self,**openparams):
        yield self.data

    def locate(self,*ordinates):
        return tuple([ax.locate(x) for x,ax in zip(ordinates,self.axes)])

    def __getitem__(self,index):
        with self.open(mode='r') as data:
            try:
                return data[index]
            except ValueError as e:
                raise IndexError(e)

    def __setitem__(self,index,value):
        with self.open() as data:
            data[index] = value

    def __iadd__(self,value):
        with self.open() as data:
            data[index] += value
    
    def __isub__(self,value):
        with self.open() as data:
            data[index] -= value
            
    def __imul__(self,value):
        with self.open() as data:
            data[index] *= value
    
    def __idiv__(self,value):
        with self.open() as data:
            data[index] /= value

    @property
    def dtype(self):
        with self.open(mode='r') as data:
            return data.dtype
    
    @property
    def values(self):
        with self.open(mode='r') as data:
            return data
            
    def interpolate(self,*axes,**kwargs):
        ndim = self.ndim
        if len(axes)!=ndim:
            raise ValueError('Expected {} dimensions'.format(ndim))
        axes = [axold.interpolate(axnew) for axold,axnew in zip(self.axes,axes)]
        axold,axnew = zip(*axes)
        return interpolate_regular(self,axold,axnew,**kwargs)
        
class NXSignalRegularGrid(RegularGrid):
    
    def __init__(self,signal,stackdim=None):
        nxdata = signal.parent
        axes = [axis.factory(values,name=name,title=attrs['title'],type='quantitative')
                for name,values,attrs in nxdata.axes]
        self.signal = signal
        super(NXSignalRegularGrid,self).__init__(axes,None,stackdim=stackdim)
    
    @contextlib.contextmanager
    def open(self,**openparams):
        with self.signal.open(**openparams) as dset:
            yield dset
            
    @property
    def values(self):
        ret = np.empty(self.shape,dtype=self.dtype)
        with self.open(mode='r') as data:
            data.read_direct(ret)
        return ret
        
class NXRegularGrid(RegularGrid):
    
    def __init__(self,nxgroup):
        # Axes in memory, data as nxfs paths
        axes = []
        signals = []
        
        if nxgroup.nxclass=='NXdata':
            it = [nxgroup]
        elif nxgroup.nxclass=='NXprocess':
            progname = nxgroup['program'].read()
            if progname!=nxfs.PROGRAM_NAME:
                raise ValueError('NXprocess from program "{}" is not known'.format(progname))
            it = nxgroup.results.iter_is_nxclass('NXdata')
        else:
            raise ValueError('{} should be an NXdata or NXprocess group'.format(nxgroup))

        for nxdata in it:
            if nxdata.islink:
                continue
            
            if not axes:
                axes = [axis.factory(values,name=name,title=attrs['title'],type='quantitative')
                        for name,values,attrs in nxdata.axes]

            for signal in nxdata.signals:
                signals.append(signal)
        
        stackdim = 0
        axnew = axis.factory(signals,type='nominal')
        axes.insert(stackdim,axnew)
        
        self.nxgroup = nxgroup
        self.signals = signals
        super(NXRegularGrid,self).__init__(axes,None,stackdim=stackdim)

    @contextlib.contextmanager
    def open(self,**openparams):
        with self.nxgroup.parent.open(**openparams) as entry:
            yield entry

    @contextlib.contextmanager
    def open_signals(self,**openparams):
        with ExitStack() as stack:
            yield [stack.enter_context(signal.open(**openparams)) for signal in self.signals]

    @property
    def signal_names(self):
        return [sig.name for sig in self.signals]

    def __getitem__(self,index):
        with self.open(mode='r') as entry:
            generator = lambda signal: entry[self.nxgroup.relpath(signal.path)]
            return indexing.getitem(generator,self.signals,index,self.ndim,shapefull=self.shape,axis=0)

    def __setitem__(self,index,value):
        with self.open(mode='r') as entry:
            selector = lambda signal: entry[self.nxgroup.relpath(signal.path)]
            indexing.setitem(selector,self.signals,index,self.ndim,value,shapefull=self.shape,axis=0,method='set')
    
    def __iadd__(self,value):
        with self.open(mode='r') as entry:
            selector = lambda signal: entry[self.nxgroup.relpath(signal.path)]
            indexing.setitem(selector,self.signals,(Ellipsis,),self.ndim,value,shapefull=self.shape,axis=0,method='add')
    
    def __isub__(self,value):
        with self.open(mode='r') as entry:
            selector = lambda signal: entry[self.nxgroup.relpath(signal.path)]
            indexing.setitem(selector,self.signals,(Ellipsis,),self.ndim,value,shapefull=self.shape,axis=0,method='sub')
            
    def __imul__(self,value):
        with self.open(mode='r') as entry:
            selector = lambda signal: entry[self.nxgroup.relpath(signal.path)]
            indexing.setitem(selector,self.signals,(Ellipsis,),self.ndim,value,shapefull=self.shape,axis=0,method='mul')
    
    def __idiv__(self,value):
        with self.open(mode='r') as entry:
            selector = lambda signal: entry[self.nxgroup.relpath(signal.path)]
            indexing.setitem(selector,self.signals,(Ellipsis,),self.ndim,value,shapefull=self.shape,axis=0,method='div')
            
    @property
    def dtype(self):
        with self.signals[0].open(mode='r') as dset:
            return dset.dtype

    @property
    def values(self):
        ret = np.empty(self.shape,dtype=self.dtype)
        for i,signal in enumerate(self.signals):
            with signal.open(mode='r') as dset:
                dset.read_direct(ret,dest_sel=(i,Ellipsis))
        return ret


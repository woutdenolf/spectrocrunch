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

from .axis import factory
from ..utils import indexing
from ..io import nxfs

class RegularGrid(with_metaclass(ABCMeta)):
    
    def __init__(self,axes):
        self.axes = axes

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
    
    @abstractmethod
    def __getitem__(self,index):
        pass

    @property
    @abstractmethod
    def values(self):
        pass
    
    @property
    @abstractmethod
    def dtype(self):
        pass
    
    def regrid(self):
        pass
    
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
            indexgen = lambda i:i,Ellipsis
            shape = self.shape[1:]
            
        elif slicedim==maxdim:
            indexgen = lambda i:Ellipsis,i
            shape = self.shape[:-1]
        else:
            a = (slice(None),)*slicedim
            b = (slice(None),)*(maxdim-slicedim)
            indexgen = lambda i:a+(i,)+b
            shape = self.shape[:slicedim]+self.shape[slicedim+1:]
        
        return shape,indexgen
    
    @contextlib.contextmanager
    def open(self,**openparams):
        yield None

class NXSignalRegularGrid(RegularGrid):
    
    def __init__(self,nxdata,signal):
        axes = [factory(values,name=name,title=attrs['title'],type='quantitative')
                for name,values,attrs in nxdata.axes]
        self.signal = signal
        super(NXSignalRegularGrid,self).__init__(axes)
    
    @contextlib.contextmanager
    def open(self,**openparams):
        with self.signal.open(**openparams) as dset:
            yield dset
    
    def __getitem__(self,index):
        with self.open(mode='r') as dset:
            try:
                return dset[index]
            except ValueError as e:
                raise IndexError(e)

    @property
    def dtype(self):
        with self.open(mode='r') as dset:
            return dset.dtype

    @property
    def values(self):
        ret = np.empty(self.shape,dtype=self.dtype)
        with self.open(mode='r') as dset:
            dset.read_direct(ret)
        return ret
        
class NXRegularGrid(RegularGrid):
    
    def __init__(self,nxgroup):
        # Axes in memory, data as nxfs paths
        axes = []
        axnew = []
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
            if not axes:
                axes = [factory(values,name=name,title=attrs['title'],type='quantitative')
                        for name,values,attrs in nxdata.axes]

            for signal in nxdata.signals:
                axnew.append((nxdata.name,signal.name))
                signals.append(signal)
                
        axnew = factory(pd.MultiIndex.from_tuples(axnew, names=['group','dataset']),
                        name='datasets',type='nominal')
        axes.insert(0,axnew)
        
        self.nxgroup = nxgroup
        self.signals = signals
        super(NXRegularGrid,self).__init__(axes)

    @contextlib.contextmanager
    def open(self,**openparams):
        with self.nxgroup.parent.open(**openparams) as entry:
            yield entry

    def __getitem__(self,index):
        with self.open(mode='r') as entry:
            generator = lambda signal: entry[self.nxgroup.relpath(signal.path)]
            return indexing.slicedstack(generator,self.signals,index,self.ndim,shapefull=self.shape,axis=0)

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


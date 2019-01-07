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

from contextlib import contextmanager
import numpy as np
try:
    from contextlib import ExitStack
except ImportError:
    from contextlib2 import ExitStack

from .regulargrid import RegularGrid
from . import axis
from . import nxresult
from ..utils import indexing
from ..io import nxfs


class NXSignalRegularGrid(RegularGrid):
    
    def __init__(self,signal,stackdim=None):
        nxdata = signal.parent
        axes = [axis.factory(values,name=name,title=attrs['title'],type='quantitative')
                for name,values,attrs in nxdata.axes]
        self.signal = signal
        super(NXSignalRegularGrid,self).__init__(axes,None,stackdim=stackdim)
    
    @contextmanager
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
        groups,axes = nxresult.regulargriddata(nxgroup)
        self.signals = [signal for group in groups.values() for signal in group]
        stackdim = 0
        axnew = axis.factory(self.signals,type='nominal')
        axes.insert(stackdim,axnew)
        self.nxentry = nxgroup.nxentry()
        super(NXRegularGrid,self).__init__(axes,None,stackdim=stackdim)

    @contextmanager
    def open(self,**openparams):
        with self.nxentry.open(**openparams) as group:
            yield group

    @contextmanager
    def open_signals(self,**openparams):
        with ExitStack() as stack:
            yield [stack.enter_context(signal.open(**openparams)) for signal in self.signals]

    @property
    def signal_names(self):
        return [sig.name for sig in self.signals]

    def __getitem__(self,index):
        with self.open(mode='r') as group:
            generator = lambda signal: group[self.nxentry.relpath(signal.path)]
            return indexing.getitem(generator,self.signals,index,self.ndim,
                    shapefull=self.shape,axis=self.stackdim)

    def __setitem__(self,index,value):
        with self.open(mode='a') as group:
            selector = lambda signal: group[self.nxentry.relpath(signal.path)]
            indexing.setitem(selector,self.signals,index,self.ndim,value,
                    shapefull=self.shape,axis=self.stackdim,method='set')
    
    def __iadd__(self,value):
        with self.open() as group:
            selector = lambda signal: group[self.nxentry.relpath(signal.path)]
            indexing.setitem(selector,self.signals,(Ellipsis,),self.ndim,value,
                    shapefull=self.shape,axis=self.stackdim,method='add')
    
    def __isub__(self,value):
        with self.open() as group:
            selector = lambda signal: group[self.nxentry.relpath(signal.path)]
            indexing.setitem(selector,self.signals,(Ellipsis,),self.ndim,value,
                    shapefull=self.shape,axis=self.stackdim,method='sub')
            
    def __imul__(self,value):
        with self.open() as group:
            selector = lambda signal: group[self.nxentry.relpath(signal.path)]
            indexing.setitem(selector,self.signals,(Ellipsis,),self.ndim,value,
                    shapefull=self.shape,axis=self.stackdim,method='mul')
    
    def __idiv__(self,value):
        with self.open() as group:
            selector = lambda signal: group[self.nxentry.relpath(signal.path)]
            indexing.setitem(selector,self.signals,(Ellipsis,),self.ndim,value,
                    shapefull=self.shape,axis=self.stackdim,method='div')
            
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

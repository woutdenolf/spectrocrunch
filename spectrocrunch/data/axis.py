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
import numpy as np
from ..utils import units
from ..utils import instance
from ..utils import listtools

class Axis(object):

    def __init__(self,params,type=None,name=None,title=None,precision=None):
        if type is None:
            type = 'quantitative'
            #quantitative: number
            #nominal: unordered categorical
            #ordinal: ordered categorical
            #temporal: date/time
            
        self.type = type
        if name is None:
            name = title
        if title is None:
            title = name
        self.name = name
        self.title = title
        self.values = params
        self.precision = precision

    @property
    def initargs(self):
        if isinstance(self._params,tuple):
            return self._params
        else:
            return self._params,
            
    @property
    def initkwargs(self):
        return {'type':self.type,'name':self.name,
                'title':self.title,'precision':self.precision}
        
    @property
    def units(self):
        return self.values.units

    @units.setter
    def units(self,value):
        self.values.ito(value)
        
    @property
    def magnitude(self):
        return self.values.magnitude
    
    def umagnitude(self,u):
        return self.values.to(u).magnitude
    
    @property
    def values(self):
        return self._values

    @values.setter
    def values(self,params):
        self._params = params
        if self.type == 'quantitative':
            params = units.asqarray(params)
        self._values = params

    @property
    def start(self):
        return self.values[0]
    
    @property
    def end(self):
        return self.values[-1]
    
    @property
    def size(self):
        return self.values.size
    
    @property
    def nsteps(self):
        return self.size-1
        
    @property
    def precision(self):
        return self._precision.to(self.units)
    
    @precision.setter
    def precision(self,value):
        if value is None:
            value = 0
        self._precision = units.Quantity(value,units=self.units)
    
    def __getitem__(self,index):
        return self.values[index]
    
    def __len__(self):
        return self.size

    def __str__(self):
        if self.units:
            return "{} ({:~})".format(self.title,self.units)
        else:
            return "{}".format(self.title)

    def __repr__(self):
        return "{}({})".format(self.name,len(self))

    def __eq__(self,other):
        if self.size!=other.size:
            return False
        diff = max(abs(self.magnitude-other.umagnitude(self.units)))
        threshold = max(self.precision,other.precision.to(self.units)).magnitude
        return diff<=self.precision.magnitude

    def __ne__(self,other):
        return not self.__eq__(other)

    def locate(self,values):
        values = units.umagnitude(values)
        x = self.magnitude
        try:
            return np.array([np.argmin(np.abs(x-v)) for v in iter(values)])          
        except TypeError:
            return np.argmin(np.abs(x-values))

    def to(self,u):
        ret = self.__class__(*self.initargs,**self.initkwargs)
        ret.units = u
        return ret
    
    def simplify(self):
        ax = self
        if self.type=='quantitative':
            kwargs = {}
            if self.size==1:
                return AxisNumber(self.start,**self.initkwargs) 
            elif self.size==2:
                return AxisRegular(self.start,self.end,1,**self.initkwargs) 
            else:
                diff = abs(np.diff(np.diff(self.magnitude)))
                if max(diff) <= self.precision.magnitude:
                    return AxisRegular(self.start,self.end,self.size-1,**self.initkwargs)
                else:
                    n = len(self)
                    ind = np.array([0] + (np.where(diff>self.precision.magnitude)[0]+1).tolist() + [n-1])
                    nseg = len(ind)-1
                    if nseg<n/3.:
                        # Limit the number of segments to 1/3 of the number of values
                        limits = self.values[ind]
                        nsteps = np.diff(ind)
                        return AxisSegments(limits,nsteps,**self.initkwargs)
                    
        return self

class _AxisRegular(Axis):

    def __init__(self,*params,**kwargs):
        kwargs['type'] = 'quantitative'
        super(_AxisRegular,self).__init__(params,**kwargs)

    @property
    def start(self):
        return self._start
    
    @start.setter
    def start(self,value):
        self._start = value
        
    @property
    def end(self):
        return self._end
        
    @end.setter
    def end(self,value):
        self._end = value
        
    @property
    def stepsize(self):
        return self._stepsize

    @stepsize.setter
    def stepsize(self,value):
        self._stepsize = value
        
    @property
    def size(self):
        return self._size

    @size.setter
    def size(self,value):
        self._size = value

    @property
    def nsteps(self):
        return self._size-1

    @nsteps.setter
    def nsteps(self,value):
        self._size = value+1
        
    @property
    def units(self):
        return self.start.units

    @units.setter
    def units(self,value):
        self.values.ito(value)
        self.start.ito(value)
        self.end.ito(value)
        self.stepsize.ito(value)

    def __repr__(self):
        return "{}(start={:~},end={:~},stepsize={:~},n={})".format(self.name,self.start,self.end,self.stepsize,len(self))
        
class AxisRegular(_AxisRegular):

    @_AxisRegular.start.setter
    def start(self,value):
        self.values = value,self.end,self.size
        
    @_AxisRegular.end.setter
    def end(self,value):
        self.values = self.start,value,self.size

    @_AxisRegular.size.setter
    def size(self,value):
        self.values = self.start,self.end,value-1
    
    @_AxisRegular.nsteps.setter
    def nsteps(self,value):
        self.values = self.start,self.end,value
        
    @_AxisRegular.values.setter
    def values(self,params):
        self._params = params
        self._start = units.Quantity(params[0])
        u = self._start.units
        self._end = units.Quantity(params[1],units=units).to(u)
        self._size = params[2]+1
        if self.size==1:
            self._stepsize = 0
        else:
            self._stepsize = (self.end-self.start)/self.nsteps
        self._values = units.Quantity(np.linspace(self.start.magnitude,self.end.magnitude,self.size),units=u)
        
class AxisRegularInc(_AxisRegular):

    @_AxisRegular.start.setter
    def start(self,value):
        self.values = value,self.stepsize,self.size

    @_AxisRegular.stepsize.setter
    def stepsize(self,value):
        self.values = self.start,value,self.size

    @_AxisRegular.size.setter
    def size(self,value):
        self.values = self.start,self.stepsize,value-1
    
    @_AxisRegular.nsteps.setter
    def nsteps(self,value):
        self.values = self.start,self.end,value
        
    @_AxisRegular.values.setter
    def values(self,params):
        self._params = params
        self._start = units.Quantity(params[0])
        u = self._start.units
        self._stepsize = units.Quantity(params[1],units=u).to(u)
        self._size = params[2]+1
        self._end = self.start + self.stepsize*self.nsteps
        self._values = np.arange(self.size)*self.stepsize+self.start

class AxisNumber(_AxisRegular):

    @_AxisRegular.start.setter
    def start(self,value):
        self.values = value,

    @_AxisRegular.values.setter
    def values(self,params):
        self._params = params
        self._start = units.Quantity(params[0])
        u = self._start.units
        self._end = self.start
        self._size = 1
        self._stepsize = units.Quantity(0,units=u)
        self._values = units.Quantity([self.start.magnitude],units=u)

    def __repr__(self):
        return "{}(={:~})".format(self.name,self.start)

class AxisSegments(Axis):
    
    def __init__(self,*params,**kwargs):
        kwargs['type'] = 'quantitative'
        super(AxisSegments,self).__init__(params,**kwargs)
        
    @property
    def start(self):
        return self.limits[0]
    
    @property
    def end(self):
        return self.limits[-1]
        
    @property
    def limits(self):
        return self._limits
    
    @limits.setter
    def limits(self,values):
        self.values = values,self.nsteps
    
    @property
    def stepsizes(self):
        return self._stepsizes
    
    @stepsizes.setter
    def stepsizes(self,values):
        limits = self.limits
        
        u = units.Quantity(values[0]).units
        if u=='dimensionless':
            u = limits.units
        func = lambda x: units.Quantity(x,units=u)
        
        nsteps = []
        for i,stepsize in enumerate(values):
            stepsize = func(stepsize)
            a,b = limits[i],limits[i+1]
            nsteps.append(int(round((b-a)/stepsize)))

        self.values = limits,nsteps

    @property
    def nsteps(self):
        return self._nsteps

    @nsteps.setter
    def nsteps(self,values):
        self.values = self.limits,values

    @staticmethod
    def mergeargs(limits,nsteps):
        nsteps = np.asarray(nsteps).tolist()+[None]
        return tuple(list(listtools.flatten(zip(limits,nsteps)))[:-1])
    
    @property
    def nsegments(self):
        return len(self._nsteps)
        
    @Axis.values.setter
    def values(self,params):
        limits,nsteps = params
        if len(limits)!=len(nsteps)+1:
            raise ValueError("Number of segments does not match the number of limits")
        self._params = params

        u = units.Quantity(limits[0]).units
        func = lambda x: units.Quantity(x,units=u).to(u).magnitude

        lmts = []
        stepsizes = []
        values = []
        for i in range(len(nsteps)):
            start,end = func(limits[i]),func(limits[i+1])
            if start >= end:
                raise ValueError("Limits must be strictly increasing")
            nspt = func(nsteps[i])
            inc = (end-start)/float(nspt)
            values += (start+inc*np.arange(nspt)).tolist()
            lmts.append(start)
            stepsizes.append(inc)
        lmts.append(end)
        values.append(end)
        
        self._stepsizes = units.Quantity(stepsizes,units=u)
        self._limits = units.Quantity(lmts,units=u)
        self._nsteps = np.asarray(nsteps)
        self._values = units.Quantity(values,units=u)
    
    def __repr__(self):
        s = ''.join(list('{:~}--({}x{:~})--'.format(lim,n,step) 
                    for lim,step,n in zip(self.limits,self.stepsizes,self.nsteps)))
        s = '{}{:~}'.format(s,self.limits[-1])
        return "{}({})".format(self.name,s)
    
def factory(data,**kwargs):
    if not isinstance(data,Axis):
        data = Axis(data,**kwargs)      
    return data.simplify()

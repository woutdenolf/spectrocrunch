# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
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

import ast
import collections
import numbers
import numpy as np
import uncertainties.core
from ..patch.pint import ureg

def isstring(x):
    return isinstance(x, basestring)

def isboollist(lst):
    try:
        return all(isinstance(i,bool) for i in lst) and len(lst)>0
    except:
        return False
        
def isarray(x):
    return isinstance(x, (list, set, frozenset, tuple, np.ndarray))

def isarray0(x):
    """Check for numpy 0-d array
    """
    if isarray(x):
        if isinstance(x,np.ndarray):
            return x.ndim==0
    return False

def isarraynot0(x):
    if isarray(x):
        if isinstance(x,np.ndarray):
            return x.ndim!=0
        else:
            return True
    return False

def isnumber(x):
    return isinstance(x, numbers.Number)
    
def isinteger(x):
    return isinstance(x, numbers.Integral)

def isscalar(x):
    return np.isscalar(x)

def isiterable(x):
    return isinstance(x, collections.Iterable)

def iscallable(x):
    return isinstance(x, collections.Callable)
    
def isquantity(x):
    return isinstance(x, ureg.Quantity)
    
def israndomvariable(x):
    # do not use asscalar!!!
    if isarray(x):
        return any(israndomvariable(z) for z in x)
    else:
        return isinstance(x,(uncertainties.core.Variable,uncertainties.core.AffineScalarFunc))

def asscalar(x):
    try:
        x = np.asscalar(x)
    except:
        pass
    return x

class _toarray(object):
    restore = {"array": lambda x:x,\
               "scalar": lambda x:x[0],\
               "array0": lambda x:np.array(x[0])}
               
    def __call__(self,x):
        if isarray(x):
            # Special case: numpy 0-d array
            if isinstance(x,np.ndarray):
                if x.ndim==0:
                    return x[np.newaxis],self.restore["array0"]
            # Create number array (possibly objects needed)
            try:
                x = np.asarray(x)
            except ValueError:
                x = np.asarray(x,dtype=object)
            return x,self.restore["array"]
        elif isquantity(x):
            u = x.units
            x,frestore = self(x.magnitude)
            func = lambda y: ureg.Quantity(y,u)
            x = np.vectorize(func,otypes=[object])(x)
            func = lambda y: ureg.Quantity(frestore(y),u)
            return x,func
        elif isnumber(x):
            return np.asarray([x]),self.restore["scalar"]
        else:
            return np.asarray([x],dtype=object),self.restore["scalar"]

asarrayf = _toarray()

def asarray(x):
    return asarrayf(x)[0]

def aslist(x):
    return asarray(x).tolist()

def asnumber(x):
    if isnumber(x):
        return x
    try:
        x = ast.literal_eval(x)
    except:
        pass
    if isnumber(x):
        return x
    else:
        try:
            return float(x)
        except:
            return np.nan
                
def arrayit(x):
    if isarray(x):
        return x
    else:
        return [x]
        

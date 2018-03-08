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

import collections
import numbers
import numpy as np
from six import string_types
import uncertainties.core
from .. import ureg

def isstring(x):
    return isinstance(x, string_types)

def isboollist(lst):
    try:
        return all(isinstance(i,bool) for i in lst) and len(lst)>0
    except:
        return False
        
def isarray(x):
    return isinstance(x, (list, set, frozenset, tuple, np.ndarray))

def isnumber(x):
    return isinstance(x, numbers.Number)
    
def isinteger(x):
    return isinstance(x, numbers.Integral)

def isscalar(x):
    return np.isscalar(x)

def isiterable(x):
    return isinstance(x, collections.Iterable)

def isquantity(x):
    return isinstance(x, ureg.Quantity)
    
def israndomvariable(x):
    if isarray(x):
        return isinstance(x.flat[0],uncertainties.core.Variable)
    else:
        return isinstance(x,uncertainties.core.Variable)

def _asarray(x,**kwargs):
    if isquantity(x):
        m = x.magnitude
        if isarray(m):
            try:
                scalar = m.ndim == 0
            except AttributeError:
                scalar = False
            if scalar:
                x = [ureg.Quantity(np.asscalar(m),x.units)]
            else:
                x = [y for y in x]
        else:
            x = [x]
        x = np.asarray(x,dtype=object,**kwargs)
    else:
        try:
            x = np.asarray(x,**kwargs)
        except ValueError:
            x = np.asarray(x,dtype=object,**kwargs)

    return x

def asarray(x,**kwargs):
    x = _asarray(x,**kwargs)
    if x.ndim == 0:
        return x[np.newaxis]
    else:
        return x
        
def asarrayf(x,**kwargs):
    x = _asarray(x,**kwargs)
    
    scalar = x.ndim == 0
    if scalar:
        # Convert to 1D array
        x = x[np.newaxis]
        func = lambda x:x[0] # not np.asscalar!!!!
    else:
        func = np.asarray
    
    return x,func

def asarrayb(x,**kwargs):
    x = _asarray(x,**kwargs)
    
    scalar = x.ndim == 0
    if scalar:
        # Convert to 1D array
        x = x[np.newaxis]
        
    return x,not scalar

def aslist(x):
    return asarray(x).tolist()

def asscalar(x):
    try:
        x = np.asscalar(x)
    except:
        pass
    return x
    



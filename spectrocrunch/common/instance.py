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

def isstring(x):
    try:
        return isinstance(x, (str, unicode))
    except:
        return isinstance(x, str)

def isboollist(lst):
    try:
        return all(isinstance(i,bool) for i in lst) and len(lst)>0
    except:
        return False
        
def isarray(x):
    return isinstance(x, (list, set, tuple, np.ndarray))

def isnumber(x):
    return isinstance(x, numbers.Number)
    
def isinteger(x):
    return isinstance(x, numbers.Integral)

def isscalar(x):
    return np.isscalar(x)

def isiterable(x):
    return hasattr(x,"__iter__")
    
def asarrayf(x,**kwargs):
    x = np.asarray(x,**kwargs)
    scalar = x.ndim == 0
    if scalar:
        # Convert to 1D array
        x = x[np.newaxis]
        func = np.asscalar
    else:
        func = np.asarray
    
    return x,func

def asarrayb(x,**kwargs):
    x = np.asarray(x,**kwargs)
    scalar = x.ndim == 0
    if scalar:
        # Convert to 1D array
        x = x[np.newaxis]
    
    return x,not scalar
    
def asarray(x,**kwargs):
    x = np.asarray(x,**kwargs)
    if x.ndim == 0:
        return x[np.newaxis]
    else:
        return x

def aslist(x):
    return asarray(x).tolist()


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

import collections
from six import string_types
import numpy as np

def _isiterable(x):
    return isinstance(x,collections.Iterable) and not isinstance(x, string_types)

def _hash_ndarray(x,numpylarge=False):
    if x.ndim==0:
        return hash(tuple(x[np.newaxis]))
    elif x.ndim==1:
        return hash(tuple(x))
    else:
        if numpylarge:
            keep,x.flags.writeable = x.flags.writeable,False
            ret = hash(x.data)
            x.flags.writeable = keep
            return ret
        else:
            return calchash(tuple(x),numpylarge=numpylarge)
            
def calchash(x,numpylarge=False):
    # Convert to iterable:
    if isinstance(x,collections.Set):
        x = sorted(x)
    elif isinstance(x,collections.MutableMapping):
        if isinstance(x,collections.OrderedDict):
            x = x.items()
        else:
            x = sorted(x.items())

    # Don't try hashing when iterable has any iterable
    xisiterable = _isiterable(x)
    tryhash = True
    if xisiterable:
        if isinstance(x,np.ndarray):
            return _hash_ndarray(x,numpylarge=numpylarge)
        elif any(_isiterable(xi) for xi in x):
            tryhash = False
        else:
            x = tuple(x)
    
    if tryhash:
        try:
            # This is only allowed to fail when iterable:
            return hash(x)
        except TypeError as e:
            if not xisiterable:
                raise e

    # Hash of the hashes of an iterable
    return hash(tuple(calchash(xi,numpylarge=numpylarge) for xi in x))
        
def hashequal(a,b,**kwargs):
    return calchash(a,**kwargs)==calchash(b,**kwargs)
    

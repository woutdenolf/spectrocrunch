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

from __future__ import division

import numpy as np
import math
import fractions
import functools
from ..utils import instance

def logscale(img):
    ret = -np.log(img/np.nanmax(img))
    ret /= np.nanmax(ret)
    return 1-ret

def sig_to_ndigits(x, sig):
    return sig-int(math.floor(math.log10(abs(x))))-1

round_ndigits = round

def ceil_ndigits(x, n):
    m = 10**n
    return math.ceil(x*m)/m

def floor_ndigits(x, n):
    m = 10**n
    return math.floor(x*m)/m
    
def round_sig(x, sig):
    return round_ndigits(x, sig_to_ndigits(x, sig))

def ceil_sig(x, sig):
    return ceil_ndigits(x, sig_to_ndigits(x, sig))

def floor_sig(x, sig):
    return floor_ndigits(x, sig_to_ndigits(x, sig))

def floatformat(x, sig):
    n = max(sig_to_ndigits(x, sig),0)
    y = "{}".format(x).split('.')
    if len(y)==2:
        n = min(n,len(y[-1]))
    return ":.0{:d}f".format(n)

def roundlist(x,max_denominator=1000000):
    x = [fractions.Fraction(a).limit_denominator(max_denominator=max_denominator) for a in x]
    denoms = [a.denominator for a in x]
    m = functools.reduce(lambda a,b: a*b//fractions.gcd(a,b), denoms)
    return np.array([int(a*m) for a in x])

def weightedsum(values,weights=None):
    if not instance.isarray(values):
        return values
    elif len(values)==1:
        return values[0]
    elif weights is None or not instance.isarray(weights):
        return np.mean(values)

    return values*weights/sum(weights)
        

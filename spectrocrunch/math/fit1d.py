# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
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

import numpy as np
import scipy.optimize
import warnings

def gaussian(x,x0,sx,A):
    return A/(np.sqrt(2*np.pi)*sx)*np.exp(-(x-x0)**2/(2*sx**2))

def errorf_gaussian(p,x,data):
    x0,sx,A = tuple(p)
    return np.ravel(gaussian(x,x0,sx,A)-data)

def guess_gaussian(x,data):
    x0i = np.argmax(data)
    x0 = x[x0i]
    sx = np.sqrt(abs((x-x0)**2*data).sum()/data.sum())
    A = data[x0]*np.sqrt(2*np.pi)*sx
    return np.array([x0,sx,A],dtype=np.float32)

def fitgaussian(x,data):
    guess = guess_gaussian(x,data)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p, success = scipy.optimize.leastsq(errorf_gaussian, guess, args=(x,data))
        success = success>0 and success<5

    return p, success


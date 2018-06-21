# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
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
"""
Statistical methods.
"""

import numpy as np

def outlierdetection(x,nsigma,pdf='normal',noutliers=0):
    
    # Detect outliers using the median-absolute-deviation
    # MAD = cte*median(x-median(x))
    # |(x-median(x))/MAD|> nsigma

    if x.size == 0:
        np.full(1,False,dtype=bool)

    # Deviation form medium
    diff = abs(x-np.median(x))
    
    # Fixed number of outliers
    if noutliers!=0:
        ret = np.full(x.size,False,dtype=bool)
        ret[(-diff).argsort()[0:noutliers]] = True
        return ret
    
    # median-absolute-deviation
    if pdf == "normal":
        MAD = 1.4826*np.median(diff)
    else:
        MAD = np.median(diff)
    
    # outliers
    if MAD == 0:
        return np.full(x.shape,False,dtype=bool)
    else:
        return diff/MAD > nsigma
    

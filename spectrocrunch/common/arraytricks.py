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
Short description of the module.
"""

import numpy as np

def tile(source,dest,reps,dirs):
    # Similar to numpy's tile, but uses stride tricks
    # dirs==0:  np.tile(source,(reps,1))
    # dirs==1:  np.tile(source,(1,reps))
    # dirs==2: diagonal tiling, not possible with np.tile

    if len(source.shape) > 2:
        raise ValueError("source dimension should be less or equal to 2.")
    if len(dest.shape) > 2:
        raise ValueError("destination dimension should be less or equal to 2.")

    if dirs == 2:
        if source.shape[0]*reps!=dest.shape[0] or source.shape[1]*reps!=dest.shape[1]:
            raise ValueError("source and destination dimensions don't correspond.")
        strides = dest.strides
        dest = np.lib.stride_tricks.as_strided(dest, shape=(reps,reps,source.shape[0],source.shape[1]), strides=(source.shape[0]*strides[0], source.shape[1]*strides[1], strides[0], strides[1]))
        ind = range(reps)
        dest[ind, ind] = source
    else:
        if dirs!=0:
            source = source.T
            dest = dest.T
        if source.shape[0]*reps!=dest.shape[0] or source.shape[1]!=dest.shape[1]:
            raise ValueError("source and destination dimensions don't correspond.")
        strides = dest.strides
        dest = np.lib.stride_tricks.as_strided(dest, shape=(reps,source.shape[0],dest.shape[1]), strides=(source.shape[0]*strides[0], strides[0], strides[1]))
        dest[:] = source

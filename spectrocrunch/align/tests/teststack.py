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

import numpy as np
import scipy.ndimage

def teststack():
    # Image and subimage sizes
    suboffset = (50,50)
    subinc = (1,2)
    nsub = 5
    maxshift = (subinc[0]*(nsub-1),subinc[1]*(nsub-1))
    subshape = (maxshift[0]*10,maxshift[1]*10)
    shape = (subshape[0] + maxshift[0] + suboffset[0]*2, subshape[1] + maxshift[1] + suboffset[1]*2)

    # Number of peaks and their width
    npeaks = 11
    nsigma = 2.

    # Large image
    n = np.round(npeaks*min(shape)/min(subshape)).astype(np.int)
    image = np.zeros(shape,dtype=np.float32)
    np.random.seed(1)
    x = shape[0]*np.random.random(n**2)
    y = shape[1]*np.random.random(n**2)
    image[x.astype(np.int), y.astype(np.int)] = 1
    image = scipy.ndimage.gaussian_filter(image, sigma=min(subshape)/(npeaks*nsigma))
    image /= np.max(image)

    # Stack of shifted subimages
    ret = np.empty(subshape+(nsub,),dtype=np.float32)
    offsets = np.empty((nsub,2),dtype=np.float32)
    for i in range(nsub):
        off = (suboffset[0] + i*subinc[0],suboffset[1] + i*subinc[1])
        offsets[i,0] = -i*subinc[0]
        offsets[i,1] = -i*subinc[1]
        ret[...,i] = image[off[0]:off[0]+subshape[0],off[1]:off[1]+subshape[1]]
    
    return ([ret,ret],offsets,2)

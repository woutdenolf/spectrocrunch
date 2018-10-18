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

def calccroproi(grid,nanval,stackdim,in_memory=False):
    """Determine crop ROI to remove slices other than stackdim which contain
    only nanval.
    
    Args:
        grid(RegularGrid):
        nanval(number):
        stackdim(int):
        in_memory(Optional(bool)):

    Returns:
        tuple or None
    """

    # Mask (True = valid pixel)
    if in_memory:
        if nanval is np.nan:
            mask = np.isnan(grid.values).sum(axis=stackdim)==0
        else:
            mask = (grid.values==nanval).sum(axis=stackdim)==0
        shape = mask.shape
    else:
        shape,indexgen = grid.sliceinfo(stackdim)
        mask = np.ones(shape,dtype=np.bool)
        for i in range(nimg):
            img = grid[indexgen(i)]
            if nanval is np.nan:
                mask &= np.logical_not(np.isnan(img))
            else:
                mask &= img != nanval

    imask = -1
    roi = []
    for igrid in range(grid.ndim):
        if igrid==stackdim:
            iroi = (0,grid.shape[stackdim])
        else:
            imask += 1
            sumdims = [i for i in range(grid.ndim-1) if i!=imask]
            indvalid = np.argwhere(mask.sum(axis=sumdims))[:,0]
            if indvalid.size==0:
                return None
            iroi = indvalid[0],indvalid[-1]+1
        roi.append(iroi)

    return roi

def convertuserroi(grid,roi,stackdim):
    """Determine crop ROI to remove slices other than stackdim which contain
    only nanval.
    
    Args:
        grid(RegularGrid):
        roi(list(2-tuple)):
        stackdim(int):

    Returns:
        list(2-tuple)
    """
    roi = list(roi)
    roi.insert(stackdim,(0,grid.shape[stackdim]))
    return roi

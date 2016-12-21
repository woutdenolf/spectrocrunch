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

import spectrocrunch.io.nexus as nexus
import math_hdf5_imagestacks_copy as m_copy

def calccroproi(stack,nanval,stackdim):
    """Determine crop ROI so that no column or row consists of only nanval

    Args:
        stack(h5py.Group): NXdata group
        nanval(number)
        stackdim(int)

    Returns:
        tuple or None
    """

    dataset = stack[stack.attrs["signal"]]

    # Stack dimensions
    if stackdim==0:
        nimg,dim1,dim2 = dataset.shape
    elif stackdim==1:
        dim1,nimg,dim2 = dataset.shape
    else:
        dim1,dim2,nimg = dataset.shape

    # Mask (True = valid pixel)
    # Do not slice by slice because of memory
    #if nanval is np.nan:
    #    mask = np.isnan(dataset[:]).sum(axis=stackdim)==0
    #else:
    #    mask = (dataset[:]==nanval).sum(axis=stackdim)==0

    mask = np.ones((dim1,dim2),dtype=np.bool)
    for i in range(nimg):
        if stackdim==0:
            img = dataset[i,...]
        elif stackdim==1:
            img = dataset[:,i,:]
        else:
            img = dataset[...,i]
        if nanval is np.nan:
            mask &= np.isnan(img)==False
        else:
            mask &= img != nanval

    # Valid row and columns (not all False)
    indvalidrow = np.argwhere(mask.sum(axis=1))
    indvalidcol = np.argwhere(mask.sum(axis=0))
    arow = indvalidrow[0]
    brow = indvalidrow[-1]+1
    acol = indvalidcol[0]
    bcol = indvalidcol[-1]+1

    # Roi to keep
    if brow-arow == dim1 and bcol-acol == dim2:
        return None
    if stackdim==0:
        roi = ((0,nimg),\
               (arow,brow),\
               (acol,bcol))
    elif stackdim==1:
        roi = ((arow,brow),\
               (0,nimg),\
               (acol,bcol))
    else:
        roi = ((arow,brow),\
               (acol,bcol),\
               (0,nimg))
    return roi

def convertuserroi(stack,roi,stackdim):

    dataset = stack[stack.attrs["signal"]]

    # Stack dimensions
    if stackdim==0:
        nimg,dim1,dim2 = dataset.shape
    elif stackdim==1:
        dim1,nimg,dim2 = dataset.shape
    else:
        dim1,dim2,nimg = dataset.shape

    # Convert user roi to crop roi
    arow = roi[0][0]
    brow = roi[0][1]
    acol = roi[1][0]
    bcol = roi[1][1]
    if arow < 0:
        arow += dim1
    if brow < 0:
        brow += dim1
    if acol < 0:
        acol += dim2
    if bcol < 0:
        bcol += dim2
    brow += 1
    bcol += 1

    # Roi to keep
    if brow-arow == dim1 and bcol-acol == dim2:
        return None
    if stackdim==0:
        roi = ((0,nimg),\
               (arow,brow),\
               (acol,bcol))
    elif stackdim==1:
        roi = ((arow,brow),\
               (0,nimg),\
               (acol,bcol))
    else:
        roi = ((arow,brow),\
               (acol,bcol),\
               (0,nimg))
    return roi

def evaluate_entire(roi,fin,varargs,retstacks):
    """ Crop image stacks based on roi
    Args:
        roi(tuple)
        fin(nexus.File)
        varargs(dict)
        fixedargs(dict)
        retstacks(list(h5py.Group))
    """
    for i in range(len(varargs)):
        grp = fin[varargs[i]]
        data = grp[grp.attrs["signal"]][roi[0][0]:roi[0][1],roi[1][0]:roi[1][1],roi[2][0]:roi[2][1]]
        dset = nexus.createNXdataSignal(retstacks[i],data=data,chunks = True)

def evaluate_sliced(roi,fin,varargs,retstacks):
    """ Crop image stacks based on roi, slice by slice
    Args:
        roi(tuple)
        fin(nexus.File)
        varargs(dict)
        fixedargs(dict)
        retstacks(list(h5py.Group))
    """
    stackdim = operation["stackdim"]

    ind = zip(np.arange(roi[stackdim][0],roi[stackdim][1]),range(roi[stackdim][1]-roi[stackdim][0]))

    for i in range(len(varargs)):
        grp = fin[varargs[i]]
        data = grp[grp.attrs["signal"]]
        dset = nexus.createNXdataSignal(retstacks[i],shape=data.shape,dtype=data.dtype,chunks = True)

        for jsource,jdest in ind:
            if stackdim==0:
                dset[jdest,...] = data[jsource,roi[1][0]:roi[1][1],roi[2][0]:roi[2][1]]
            elif stackdim==1:
                dset[:,jdest,:] = data[roi[0][0]:roi[0][1],jsource,roi[2][0]:roi[2][1]]
            else:
                dset[...,jdest] = data[roi[0][0]:roi[0][1],roi[1][0]:roi[1][1],jsource]

def evaluate(operation,fin,varargs,retstacks,axes):
    """ Crop image stacks based on a nanvalue or on manual range
    Args:
        operation(dict)
        fin(nexus.File)
        varargs(dict)
        fixedargs(dict)
        retstacks(list(h5py.Group))
    """
    axesdata = []

    reference = [s for s in varargs if s.endswith(operation["reference set"])]
    iref = varargs.index(reference[0])

    roi = None
    if "nanval" in operation:
        roi = calccroproi(fin[varargs[iref]],operation["nanval"],operation["stackdim"])
    elif "roi" in operation:
        roi = convertuserroi(fin[varargs[iref]],operation["roi"],operation["stackdim"])

    # Add (cropped) data to the NXdata group
    if roi is None:
        m_copy.evaluate(operation,fin,varargs,retstacks)
    else:
        bdone = False
        if "sliced" in operation:
            if operation["sliced"]:
                evaluate_sliced(roi,fin,varargs,retstacks)
                bdone = True
        if not bdone:
            evaluate_entire(roi,fin,varargs,retstacks)

        axesdata = [None]*len(axes)
        for i in range(len(axes)):
            axesdata[i] = fin[axes[i]["fullname"]][roi[i][0]:roi[i][1]]

    return axesdata

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

import numpy as np
from . import regulargrid
from . import nxutils

def calccroproi(grid,nanval,stackdim,sliced=True):
    """Determine crop ROI to remove slices other than stackdim which contain
    only nanval.
    
    Args:
        grid(RegularGrid):
        nanval(number):
        stackdim(int):
        sliced(Optional(bool)):

    Returns:
        tuple or None
    """

    # Mask (True = valid pixel)
    if sliced:
        shape,indexgen = grid.sliceinfo(stackdim)
        mask = np.ones(shape,dtype=np.bool)
        for i in range(grid.shape[stackdim]):
            img = grid[indexgen(i)]
            if nanval is np.nan:
                mask &= np.logical_not(np.isnan(img))
            else:
                mask &= img != nanval
    else:
        if nanval is np.nan:
            mask = np.isnan(grid.values).sum(axis=stackdim)==0
        else:
            mask = (grid.values==nanval).sum(axis=stackdim)==0
        shape = mask.shape

    imask = -1
    roi = []
    for igrid in range(grid.ndim):
        if igrid==stackdim:
            iroi = (0,grid.shape[stackdim])
        else:
            imask += 1
            sumdims = tuple([i for i in range(grid.ndim-1) if i!=imask])
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

def execute(grid,roi,parameters):
    """Apply ROI
    
    Args:
        grid(RegularGrid):
        roi(list(2-tuple)):
        parameters(dict):

    Returns:
        nxfs.Path
    """
    # New nxprocess (return when already exists)
    process,notnew = nxutils.next_process([grid.nxgroup],parameters)
    if notnew:
        return process
    results = process.results
        
    # Create new axes (if needed)
    positioners = grid.nxgroup.positioners()
    axes = []
    gaxes = grid.axes
    gaxes.pop(grid.stackdim)
    for ax,(a,b) in zip(gaxes,roi):
        if ax.size==b-a:
            axes.append(ax.name)
        else:
            name = '{}_{}'.format(ax.name,parameters['method'])
            positioners.add_axis(name,ax[a:b],title=ax.title)
            axes.append(name)


    
    

    # Cropped signals
    stackdim = parameters["stackdim"]
    shape = tuple([b-a for a,b in roi])
    indexcrop = [slice(a,b) for a,b in roi]
    indexout = [slice(None)]*len(roi)
    nslice = grid.shape[stackdim]
    dtype = grid.dtype
    sliced = parameters["sliced"]
    for signalin in grid.signals:
        nxdata = results[signalin.parent.name]
        bnew = not nxdata.exists
        if bnew:
            nxdata = results.nxdata(name=signalin.parent.name)

        with signalin.open() as dsetin:
            if sliced:
                signalout = nxdata.add_signal(signalin.name,shape=shape,
                                              dtype=dtype,chunks=True)
                with signalout.open() as dsetout:
                    for i in range(nslice):
                        indexcrop[stackdim] = i
                        indexout[stackdim] = i
                        dsetout[tuple(indexout)] = dsetin[tuple(indexcrop)]
            else:
                signalout = nxdata.add_signal(signalin.name,data=dsetin[tuple(indexcrop)],
                                              chunks=True)
    
        if bnew: 
            nxdata.set_axes(*axes)
        
    return process
    
def crop(grid,parameters):
    # Parameters with defaults
    parameters["sliced"] = parameters.get("sliced",True)
    parameters["stackdim"] = parameters.get("stackdim",0)
    
    # Signal to determine the ROI
    reference = parameters["reference"]
    ax = grid.axes[0]
    ind = np.array([s.path.endswith(reference) for s in ax])
    if ind.any():
        ind = np.where(ind)[0][-1]
        refgrid = regulargrid.NXSignalRegularGrid(ax[ind])
    else:
        raise ValueError('Reference "{}" not present in {}'.format(ref,ax))
    
    # Determine ROI
    roi = None
    if "nanval" in parameters:
        roi = calccroproi(refgrid,parameters["nanval"],parameters["stackdim"],
                          sliced=parameters["sliced"])
    elif "roi" in parameters:
        roi = convertuserroi(refgrid,parameters["roi"],parameters["stackdim"])
    else:
        raise ValueError('Specify either "nanval" or "roi"')
    
    # Apply ROI
    if roi:
        return execute(grid,roi,parameters)
    else:
        return None
        

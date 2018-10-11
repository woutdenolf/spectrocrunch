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

import json
import numpy as np
import h5py
import fabio
import glob
import re

from ..utils.dict import defaultdict
from ..io import nexus

def execrebin(img,rebin):
    """
    Args:
        img(np.array): 2 dimensions
        rebin(array-like)
    Returns:
        np.array: rebinned image
    """
    shaperebin = (data.shape[0]//rebin[0],data.shape[1]//rebin[1])
    view = data[:shaperebin[0] * rebin[0],:shaperebin[1] * rebin[1]]
    view = view.reshape(shaperebin[0], rebin[0], shaperebin[1], rebin[1])
    return view.sum(axis=3).sum(axis=1)

def execroi(img,roi):
    """
    Args:
        img(np.array): 2 dimensions
        roi(2-tuple(2-tuple))
    Returns:
        np.array: cropped image
    """
    return img[roi[0][0]:roi[0][1],roi[1][0]:roi[1][1]]

def procraw(img,config):
    roi = config["roi"]
    if roi is not None:
        img = execroi(img,roi)

    rebin = config["rebin"]
    if rebin[0] > 1 or rebin[1] > 1:
        img = execrebin(img,rebin)

    return img

def darklibrary(config):
    """
    Args:
        config(dict)
    Returns:
        spectrocrunch.utils.dict.defaultdict: dictionary of dark frames for particular exposures
    """

    darkfiles = []
    for f in config["darklist"]:
        darkfiles += glob.glob(f)

    # Labels
    frametimelabel = config["frametimelabel"]
    frametimedefault = str(config["frametimedefault"])
    dtype = eval(config["dtype"])
    nflabel = config["nbframeslabel"]

    # Library
    dark = defaultdict()
    dark.setdefaultfactory(lambda frametime: {"data":config["darkcurrentzero"] + config["darkcurrentgain"]*float(frametime),"nframes":1})
    
    for f in darkfiles:
        fh = fabio.open(f)
        h = fh.header

        # Raw data
        data = fh.data.astype(dtype)
        data = procraw(data,config)
        
        # Frame time
        if frametimelabel in h:
            frametime = h[frametimelabel]
        else:
            frametime = frametimedefault

        # Number of frames
        if nflabel in h:
            nframes = dtype(h[nflabel])
        else:
            nframes = dtype(1)
    
        # Add to library
        if frametime in dark:
            dark[frametime]["data"] += data
            dark[frametime]["nframes"] += nframes
        else:
            dark[frametime] = {"data":data,"nframes":nframes}

    return dark

def getroi(roilabel,h,fh):
    """ roi = [row0,col0,nrow,ncol] or [y0,x0,ny,nx]
    """

    # ROI from label
    roi = []
    if roilabel in h:
        roi = [int(s) for s in re.split('[<>,\-x ]+',h[roilabel]) if len(s)!=0]
        roi = [roi[1],roi[0],roi[3],roi[2]] # camera takes x as the first dimension, python takes y

    # ROI from dimensions
    if len(roi)!=4:
        roi = [0,0,int(fh.dim2),int(fh.dim1)] # edf: dim2 => rows, dim1 => cols

    return roi

def axesvalues(roi,config):
    x0 = roi[1]
    y0 = roi[0]
    nx = roi[3]
    ny = roi[2]

    # Sub region
    if config["roi"] is not None:
        ax = config["roi"][1][0]
        bx = config["roi"][1][1]
        ay = config["roi"][0][0]
        by = config["roi"][0][1]
        if ax is None:
            ax = 0
        if bx is None:
            bx = nx
        if ay is None:
            ay = 0
        if by is None:
            by = ny
        if ax < 0:
            ax += nx
        if bx < 0:
            bx += nx
        if ay < 0:
            ay += ny
        if by < 0:
            by += ny

        x0 += ax
        y0 += ay
        nx = bx-ax
        ny = by-ay
       
    # Rebin
    dx = config["rebin"][1]
    dy = config["rebin"][0]

    nx //= dx
    ny //= dy

    x0 += (dx-1)/2.
    y0 += (dy-1)/2.
    
    x1 = x0 + (nx-1)*dx
    y1 = y0 + (ny-1)*dy

    # axes in pixels
    col = np.linspace(x0,x1,nx)
    row = np.linspace(y0,y1,ny)

    return row,col

def dataflatlibrary(config):
    """Separate and sort data and flat fields

    Args:
        config(dict): description
        
    Returns:
        data(dict):
        flat1(dict):
        flat2(dict|None):
        keyindices(list): indices for sorting
        stackaxes(dict): [{"name":"name1","data":np.array},
        {"name":"name2","data":np.array},
        {"name":"name3","data":np.array}]
    """

    datafiles = []
    for f in config["datalist"]:
        datafiles += glob.glob(f)

    flatfiles = []
    for f in config["flatlist"]:
        flatfiles += glob.glob(f)

    # Labels
    stacklabel = config["stacklabel"]
    roilabel = config["roilabel"]

    # Data and flat dictionaries
    roi = None
    data = {}
    flat = {}
    for f in datafiles:
        fh = fabio.open(f)
        h = fh.header

        # Stack label
        if stacklabel not in h:
            continue
        key = h[stacklabel]

        # ROI
        tmp = getroi(roilabel,h,fh)
        if roi is None:
            roi = tmp
        else:
            if tmp != roi:
                continue

        # Add to library
        if key in data:
            data[key].append(f)
        else:
            data[key] = [f]
            flat[key] = []

    for f in flatfiles:
        fh = fabio.open(f)
        h = fh.header

         # Stack label
        if stacklabel not in h:
            continue
        key = h[stacklabel]

        # ROI
        tmp = getroi(roilabel,h,fh)
        if roi is None:
            roi = tmp
        else:
            if tmp != roi:
                continue

        # Add to library
        if key in flat:
            flat[key].append(f)

    if len(data)==0:
        raise IOError("Not data files found ({})".format(config["datalist"][0]))

    # Remove entries with missing images
    #ndata = map(lambda x:len(data[x]),data)
    #ndatafiles = max(ndata,key=ndata.count)
    #nflat = map(lambda x:len(flat[x]),flat)
    #nflatfiles = max(nflat,key=nflat.count)
    #tmp1 = {k:data[k] for k in data if len(data[k])==ndatafiles and len(flat[k])==nflatfiles}
    #tmp2 = {k:flat[k] for k in flat if len(data[k])==ndatafiles and len(flat[k])==nflatfiles}
    #data = tmp1
    #flat = tmp2
    
    # Split up flat images
    if config["beforeafter"]:
        flat1 = {}
        flat2 = {}
        for k in flat:
            files = sorted(flat[k])
            n = len(files)
            if n==1:
                flat1[k] = files
                flat2[k] = []
            else:
                flat1[k] = files[0:n//2]
                flat2[k] = files[n//2:]
    else:
        flat1 = flat
        flat2 = None

    # Axes
    stackvalues = np.array(map(lambda x:np.float32(x),data.keys()))
    keyindices = np.argsort(stackvalues)
    stackvalues = stackvalues[keyindices]
    row,col = axesvalues(roi,config)

    stackdim,imgdim = dimensions(config)
    stackaxes = [None]*3
    stackaxes[imgdim[0]] = {"name":"row","data":row}
    stackaxes[imgdim[1]] = {"name":"col","data":col}
    stackaxes[stackdim] = {"name":stacklabel,"data":stackvalues}

    # Result
    return data,flat1,flat2,keyindices,stackaxes

def dimensions(config):
    stackdim = config["stackdim"]
    if stackdim == 0:
        imgdim = [1,2]
    elif stackdim == 1:
        imgdim = [0,2]
    else:
        imgdim = [0,1]
    return stackdim,imgdim

def getsingleimage(filename,darklib,config):
    """ Get image and corresponding information (dark, expo time, nframes)
    """
    # Labels
    frametimelabel = config["frametimelabel"]
    frametimedefault = str(config["frametimedefault"])
    dtype = eval(config["dtype"])
    nflabel = config["nbframeslabel"]

    fh = fabio.open(filename)
    h = fh.header

    # Raw data
    data = fh.data.astype(dtype)
    data = procraw(data,config)

    # Frame time
    if frametimelabel in h:
        frametime = h[frametimelabel]
    else:
        frametime = frametimedefault

    # Number of frames
    if nflabel in h:
        nframes = dtype(h[nflabel])
    else:
        nframes = dtype(1)

    return data,frametime,nframes,darklib[frametime]

def getnormalizedimage(fileslist,darklib,config):
    """ Get dark subtracted images from a list of files with intensity in DU/sec
        
        img = (img1 - nf1*dark1) + (img2 - nf2*dark2) + ...
        time = nf1*tframe1 + nf2*tframe2 + ...
        img /= time
    """

    # Labels
    frametimelabel = config["frametimelabel"]
    frametimedefault = str(config["frametimedefault"])
    dtype = eval(config["dtype"])
    nflabel = config["nbframeslabel"]

    img = None
    time = None
    for f in fileslist:
        fh = fabio.open(f)
        h = fh.header

        # Raw data
        data = fh.data.astype(dtype)
        data = procraw(data,config)

        # Frame time
        if frametimelabel in h:
            frametime = h[frametimelabel]
        else:
            frametime = frametimedefault

        # Number of frames
        if nflabel in h:
            nframes = dtype(h[nflabel])
        else:
            nframes = dtype(1)

        data -= darklib[frametime]["data"]/darklib[frametime]["nframes"]*nframes
        if img is None:
            img = data
            time = dtype(frametime)*nframes
        else:
            img += data
            time += dtype(frametime)*nframes

    img /= time
    return img

def create_hdf5_imagestacks(jsonfile):
    """Convert transmission images to an HDF5 file:
        groups which contain NXdata classes
        3 axes datasets on the main level

    Returns:
        stacks: {"counters":{"name1":lstack1,"name2":lstack2,...},
                 "det0":{"name3":lstack3,"name4":lstack4,...},
                 "det1":{"name3":lstack5,"name4":lstack6,...},...}
                 lstack: an image stack given as an NXdata path

        axes: [{"name":"name1","fullname":"/axes/name1/data"},
               {"name":"name2","fullname":"/axes/name2/data"},
               {"name":"name3","fullname":"/axes/name3/data"}]
    """

    # This is the structure at ID21 (not imposed):
    # data: path/radix_zonenumber_energynumber_repeat.edf
    # flat: path/radix_zonenumber_energynumber_repeat.edf
    # path/radix_dark_exposuretime_repeat.edf

    # Processing configuration
    with open(jsonfile,'r') as f:
        config = json.load(f)

    # Dark, data and flat libraries
    darklib = darklibrary(config)
    data,flat1,flat2,keyindices,stackaxes = dataflatlibrary(config)

    # Prepare stack dimensions
    dim = [0]*3
    stackdim,imgdim = dimensions(config)

    # Open file and create detector group
    f = nexus.File(config["hdf5output"],mode='w')
    grp = nexus.newNXentry(f,"detector0")

    # Save stack axes values
    axes = nexus.createaxes(f,stackaxes)

    # Create NXdata groups for transmission and flat-field stacks
    nxdatasample = nexus.newNXdata(grp,"sample","")
    nxdataflat1 = None
    nxdataflat2 = None
    if not config["normalize"]:
        nxdataflat1 = nexus.newNXdata(grp,"flat1","")
        if flat2 is not None:
            nxdataflat2 = nexus.newNXdata(grp,"flat2","")

    # Loop over the images
    keys = data.keys()
    dtype = eval(config["dtype"])
    for i in range(len(keyindices)):
        key = keys[keyindices[i]]
    
        # Get data
        img = getnormalizedimage(data[key],darklib,config)

        # Normalize
        if config["normalize"]:
            flat = getnormalizedimage(flat1[key],darklib,config)
            if flat2 is not None:
                flat += getnormalizedimage(flat2[key],darklib,config)
                flat /= 2.
            img = -np.log(img/flat)

        # Allocate datasets
        if i==0:
            dim[imgdim[0]] = img.shape[0]
            dim[imgdim[1]] = img.shape[1]
            dim[stackdim] = len(data)
            dsetsample = nexus.createNXdataSignal(nxdatasample,shape=dim,chunks = True,dtype = dtype)
            if nxdataflat1 is not None:
                dsetflat1 = nexus.createNXdataSignal(nxdataflat1,shape=dim,chunks = True,dtype = dtype)
                if nxdataflat2 is not None:
                    dsetflat2 = nexus.createNXdataSignal(nxdataflat2,shape=dim,chunks = True,dtype = dtype)

        if stackdim == 0:
            dsetsample[i,...] = img
            if nxdataflat1 is not None:
                dsetflat1[i,...] = getnormalizedimage(flat1[key],darklib,config)
                if nxdataflat2 is not None:
                    if len(flat2[key])==0:
                        dsetflat2[i,...] = dsetflat1[i,...]
                    else:
                        dsetflat2[i,...] = getnormalizedimage(flat2[key],darklib,config)
        elif stackdim == 1:
            dsetsample[:,i,:] = img
            if nxdataflat1 is not None:
                dsetflat1[:,i,:] = getnormalizedimage(flat1[key],darklib,config)
                if nxdataflat2 is not None:
                    if len(flat2[key])==0:
                        dsetflat2[:,i,:] = dsetflat1[:,i,:]
                    else:
                        dsetflat2[:,i,:] = getnormalizedimage(flat2[key],darklib,config)
        else:
            dsetsample[...,i] = img
            if nxdataflat1 is not None:
                dsetflat1[...,i] = getnormalizedimage(flat1[key],darklib,config)
                if nxdataflat2 is not None:
                    if len(flat2[key])==0:
                        dsetflat2[...,i] = dsetflat1[...,i]
                    else:
                        dsetflat2[...,i] = getnormalizedimage(flat2[key],darklib,config)

    # Stack dict and link axes
    if nxdataflat1 is None:
        stacks = {"detector0":{"sample":nxdatasample.name}}
        nexus.linkaxes(f,axes,[nxdatasample])
    else:
        if nxdataflat2 is None:
            stacks = {"detector0":{"sample":nxdatasample.name,"flat1":nxdataflat1.name}}
            nexus.linkaxes(f,axes,[nxdatasample,nxdataflat1])
        else:
            stacks = {"detector0":{"sample":nxdatasample.name,"flat1":nxdataflat1.name,"flat2":nxdataflat2.name}}
            nexus.linkaxes(f,axes,[nxdatasample,nxdataflat1,nxdataflat2])

    # Save stackinfo
    #stackinfogrp = nexus.newNXentry(f,"stackinfo")
    #for k in stackinfo:
    #    stackinfogrp[k] = stackinfo[k]

    # Add processing info
    #nexus.addinfogroup(f,"fromraw",config)
    nexus.addinfogroup(f,"fromraw",{"config":jsonfile,"pixel unit":"DU/sec","dark current":"subtracted"})

    f.close()

    return stacks,axes

if __name__ == '__main__':
    import sys
    if len(sys.argv)>=2:
        create_hdf5_imagestacks(sys.argv[1])


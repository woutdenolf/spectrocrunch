# -*- coding: utf-8 -*-

import json
import numpy as np
import h5py
import fabio
import glob
import re

from ..utils.dict import defaultdict


def execrebin(img, rebin):
    """
    Args:
        img(np.array): 2 dimensions
        rebin(array-like)
    Returns:
        np.array: rebinned image
    """
    shaperebin = (data.shape[0]//rebin[0], data.shape[1]//rebin[1])
    view = data[:shaperebin[0] * rebin[0], :shaperebin[1] * rebin[1]]
    view = view.reshape(shaperebin[0], rebin[0], shaperebin[1], rebin[1])
    return view.sum(axis=3).sum(axis=1)


def execroi(img, roi):
    """
    Args:
        img(np.array): 2 dimensions
        roi(2-tuple(2-tuple))
    Returns:
        np.array: cropped image
    """
    return img[roi[0][0]:roi[0][1], roi[1][0]:roi[1][1]]


def procraw(img, parameters):
    roi = parameters["roi"]
    if roi is not None:
        img = execroi(img, roi)

    rebin = parameters["rebin"]
    if rebin[0] > 1 or rebin[1] > 1:
        img = execrebin(img, rebin)

    return img


def darklibrary(parameters):
    """
    Args:
        parameters(dict)
    Returns:
        spectrocrunch.utils.dict.defaultdict: dictionary of dark frames for particular exposures
    """

    darkfiles = []
    for f in parameters["darklist"]:
        darkfiles += glob.glob(f)

    # Labels
    frametimelabel = parameters["frametimelabel"]
    frametimedefault = str(parameters["frametimedefault"])
    dtype = eval(parameters["dtype"])
    nflabel = parameters["nbframeslabel"]

    # Library
    dark = defaultdict()
    dark.setdefaultfactory(lambda frametime: {
                           "data": parameters["darkcurrentzero"] + parameters["darkcurrentgain"]*float(frametime), "nframes": 1})

    for f in darkfiles:
        fh = fabio.open(f)
        h = fh.header

        # Raw data
        data = fh.data.astype(dtype)
        data = procraw(data, parameters)

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
            dark[frametime] = {"data": data, "nframes": nframes}

    return dark


def getroi(roilabel, h, fh):
    """
    Returns:
        list: [row0,col0,nrow,ncol] or [y0,x0,ny,nx]
    """

    # ROI from label
    roi = []
    if roilabel in h:
        roi = [int(s)
               for s in re.split('[<>,\-x ]+', h[roilabel]) if len(s) != 0]
        # camera takes x as the first dimension, python takes y
        roi = [roi[1], roi[0], roi[3], roi[2]]

    # ROI from dimensions
    if len(roi) != 4:
        # edf: dim2 => rows, dim1 => cols
        roi = [0, 0, int(fh.dim2), int(fh.dim1)]

    return roi


def axesvalues(roi, parameters):
    """
    Args:
        roi(list): y0,x0,ny,nx
        parameters(dict)
    Returns:
        tuple(np.ndarray)
    """
    x0 = roi[1]
    y0 = roi[0]
    nx = roi[3]
    ny = roi[2]

    # Sub region
    if parameters["roi"] is not None:
        ax = parameters["roi"][1][0]
        bx = parameters["roi"][1][1]
        ay = parameters["roi"][0][0]
        by = parameters["roi"][0][1]
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
    dx = parameters["rebin"][1]
    dy = parameters["rebin"][0]

    nx //= dx
    ny //= dy

    x0 += (dx-1)/2.
    y0 += (dy-1)/2.

    x1 = x0 + (nx-1)*dx
    y1 = y0 + (ny-1)*dy

    # axes in pixels
    col = np.linspace(x0, x1, nx)
    row = np.linspace(y0, y1, ny)

    return row, col


def dataflatlibrary(parameters):
    """Separate and sort data and flat fields

    Args:
        parameters(dict): description

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
    for f in parameters["datalist"]:
        datafiles += glob.glob(f)

    flatfiles = []
    for f in parameters["flatlist"]:
        flatfiles += glob.glob(f)

    # Labels
    stacklabel = parameters["stacklabel"]
    roilabel = parameters["roilabel"]

    # Stack values
    datastackvalues = parameters["datastackvalues"]
    if datastackvalues:
        if not len(datastackvalues) == len(datafiles):
            datastackvalues = []
    if not datastackvalues:
        datastackvalues = [None] * len(datafiles)
    flatstackvalues = parameters["flatstackvalues"]
    if flatstackvalues:
        if not len(flatstackvalues) == len(flatfiles):
            flatstackvalues = []
    if not flatstackvalues:
        flatstackvalues = [None] * len(flatfiles)

    # Data and flat dictionaries
    roi = None
    data = {}
    flat = {}
    for f, stackvalue in zip(datafiles, datastackvalues):
        fh = fabio.open(f)
        h = fh.header

        # Stack label
        if stackvalue is None:
            if stacklabel not in h:
                continue
            key = h[stacklabel]
        else:
            key = stackvalue

        # ROI
        tmp = getroi(roilabel, h, fh)
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

    for f, stackvalue in zip(flatfiles, flatstackvalues):
        fh = fabio.open(f)
        h = fh.header

        # Stack label
        if stackvalue is None:
            if stacklabel not in h:
                continue
            key = h[stacklabel]
        else:
            key = stackvalue

        # ROI
        tmp = getroi(roilabel, h, fh)
        if roi is None:
            roi = tmp
        else:
            if tmp != roi:
                continue

        # Add to library
        if key in flat:
            flat[key].append(f)

    if len(data) == 0:
        raise IOError("No data files found ({})".format(
            parameters["datalist"][0]))

    # Remove entries with missing images
    #ndata = list(map(lambda x:len(data[x]),data))
    #ndatafiles = max(ndata,key=ndata.count)
    #nflat = list(map(lambda x:len(flat[x]),flat))
    #nflatfiles = max(nflat,key=nflat.count)
    #tmp1 = {k:data[k] for k in data if len(data[k])==ndatafiles and len(flat[k])==nflatfiles}
    #tmp2 = {k:flat[k] for k in flat if len(data[k])==ndatafiles and len(flat[k])==nflatfiles}
    #data = tmp1
    #flat = tmp2

    # Split up flat images
    if parameters["flatbeforeafter"]:
        flat1 = {}
        flat2 = {}
        for k in flat:
            files = sorted(flat[k])
            n = len(files)
            if n == 1:
                flat1[k] = files
                flat2[k] = []
            else:
                flat1[k] = files[0:n//2]
                flat2[k] = files[n//2:]
    else:
        flat1 = flat
        flat2 = None

    # Axes
    stackvalues = np.array(list(map(lambda x: np.float32(x), data.keys())))
    keyindices = np.argsort(stackvalues)
    stackvalues = stackvalues[keyindices]
    row, col = axesvalues(roi, parameters)

    stackdim = parameters["stackdim"]
    imgdim = [i for i in range(3) if i != stackdim]
    stackaxes = [None]*3
    stackaxes[imgdim[0]] = {"name": "row", "data": row}
    stackaxes[imgdim[1]] = {"name": "col", "data": col}
    stackaxes[stackdim] = {"name": stacklabel, "data": stackvalues}

    # Result
    return data, flat1, flat2, keyindices, stackaxes


def getsingleimage(filename, darklib, parameters):
    """ Get image and corresponding information (dark, expo time, nframes)

    Returns:
        data(np.ndarray): image
        frametime(num): exposure time per frame
        nframes(num): number of frames
        dark(dict): {"data":...,"nframes":...}
    """
    # Labels
    frametimelabel = parameters["frametimelabel"]
    frametimedefault = str(parameters["frametimedefault"])
    dtype = eval(parameters["dtype"])
    nflabel = parameters["nbframeslabel"]

    fh = fabio.open(filename)
    h = fh.header

    # Raw data
    data = fh.data.astype(dtype)
    data = procraw(data, parameters)

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

    return data, frametime, nframes, darklib[frametime]


def getnormalizedimage(fileslist, darklib, parameters):
    """ Get dark subtracted images from a list of files with intensity in DU/sec

        img = (img1 - nf1*dark1) + (img2 - nf2*dark2) + ...
        time = nf1*tframe1 + nf2*tframe2 + ...
        img /= time
    """

    # Labels
    frametimelabel = parameters["frametimelabel"]
    frametimedefault = str(parameters["frametimedefault"])
    dtype = eval(parameters["dtype"])
    nflabel = parameters["nbframeslabel"]

    img = None
    time = None
    for f in fileslist:
        fh = fabio.open(f)
        h = fh.header

        # Raw data
        data = fh.data.astype(dtype)
        data = procraw(data, parameters)

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

        data -= darklib[frametime]["data"] / \
            darklib[frametime]["nframes"]*nframes
        if img is None:
            img = data
            time = dtype(frametime)*nframes
        else:
            img += data
            time += dtype(frametime)*nframes

    img /= time
    return img

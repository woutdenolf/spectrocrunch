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

import logging
import os
import json
import numpy as np

import spectrocrunch.common.timing as timing
from spectrocrunch.fullfield.create_hdf5_imagestacks import create_hdf5_imagestacks as makestacks
from spectrocrunch.h5stacks.get_hdf5_imagestacks import get_hdf5_imagestacks as getstacks

import spectrocrunch.io.nexus as nexus

from .proc_normalize import execute as normalize
from .proc_align import execute as align
from .proc_crop import execute as execcrop
from . proc_common import defaultstack
from . proc_common import flattenstacks

def createconfig_pre(sourcepath,destpath,radix,ext,rebin,roi,stackdim):

    if not isinstance(sourcepath,list):
        sourcepath = [sourcepath]
    if not isinstance(radix,list):
        radix = [radix]

    config = {
        # EDF header
        "frametimelabel" : "exposure_time",
        "frametimedefault" : 0.,
        "roilabel" : "roi",
        "nbframeslabel" : "nb_frames",
        "stacklabel" : "energy",

        # Defaults
        "dtype": "np.float32", # must be floating point!
        "darkcurrentzero" : 90.,
        "darkcurrentgain" : 1.,
        "stackdim" : stackdim,

        # Data
        "darklist" : map(lambda xy: os.path.join(xy[0],xy[1]+"*_dark_*.edf"),zip(sourcepath,radix)),
        "datalist" : map(lambda xy: os.path.join(xy[0],xy[1]+"*_data_*.edf"),zip(sourcepath,radix)),
        "flatlist" : map(lambda xy: os.path.join(xy[0],xy[1]+"*_ref_*.edf"),zip(sourcepath,radix)),
        "beforeafter" : True,

        # Output
        "hdf5output": os.path.join(destpath,radix[0]+ext+".h5"),
        "roi": roi,
        "rebin":rebin}

    # Create configuration file
    if not os.path.exists(destpath):
        os.makedirs(destpath)
    jsonfile = os.path.join(destpath,radix[0]+ext+".json")
    with open(jsonfile,'w') as f:
        json.dump(config,f,indent=2)

    return jsonfile,config["hdf5output"]

def process(sourcepath,destpath,radix,ext,rebin,alignmethod,\
        skippre=False,skipnormalization=False,skipalign=False,\
        roiraw=None,roialign=None,roiresult=None,\
        refimageindex=None,crop=False,plot=True):

    logger = logging.getLogger(__name__)
    T0 = timing.taketimestamp()

    stackdim = 2
    bsamefile = False

    cropalign = crop

    # Image stack
    jsonfile, h5file = createconfig_pre(sourcepath,destpath,radix,ext,rebin,roiraw,stackdim)
    if skippre:
        stacks, axes = getstacks(h5file,["detector0"])
    else:
        logger.info("Creating image stacks ...")
        stacks, axes = makestacks(jsonfile)

        #stacks2, axes2 = getstacks(h5file,["detector0"])
        #assert(axes == axes2)
        #assert(stacks == stacks2)

    # Convert stack dictionary to stack list
    stacks = flattenstacks(stacks)

    # Default group
    default = "sample"
    defaultstack(h5file,stacks,default)

    # Groups that don't change and need to be copied
    copygroups = None

    # I0 normalization
    if skipnormalization:
        file_normalized, Ifn_stacks,Ifn_axes = h5file,stacks,axes
    else:
        if any("flat2" in s for s in stacks):
            snorm = ["flat1","flat2"]
        else:
            snorm = "flat1"
        file_normalized, Ifn_stacks,Ifn_axes = normalize(h5file,stacks,axes,copygroups,bsamefile,default,snorm,snorm,stackdim=stackdim,copyskipped=False)

    # Alignment
    if alignmethod is None or skipalign:
        timing.printtimeelapsed(T0,logger)
        exit()
    else:
        alignreference = "sample"
        file_aligned, aligned_stacks, aligned_axes = align(file_normalized, Ifn_stacks, Ifn_axes, copygroups, bsamefile, default,\
            alignmethod, alignreference, refimageindex, cropalign, roialign, plot, stackdim)

    # Crop
    if roiresult is not None:
        cropinfo = {"roi":roiresult,"stackdim":stackdim,"reference set":"sample"}
        execcrop(h5file, stacks, axes, copygroups, bsamefile, default, cropinfo)

    timing.printtimeelapsed(T0,logger)



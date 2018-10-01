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

from ..common import timing
from ..fullfield.create_hdf5_imagestacks import create_hdf5_imagestacks as makestacks
from ..h5stacks.get_hdf5_imagestacks import get_hdf5_imagestacks as getstacks

from ..io import nexus

from .proc_math import execute as normalizefunc
from .proc_align import execute as align
from .proc_crop import execute as execcrop
from .proc_common import defaultstack
from .proc_common import flattenstacks

def createconfig(sourcepath,destpath,radix,**kwargs):
    rebin = kwargs.get("rebin",(1,1))
    roiraw = kwargs.get("roiraw",None)
    stackdim = kwargs.get("stackdim",2)
    skipnormalization = kwargs.get("skipnormalization",False)
    normalizeonload = kwargs.get("normalizeonload",True)
    normalize = normalizeonload and not skipnormalization
    extout = kwargs.get("extout","").replace(".","_")

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
        "beforeafter" : True, # split up flat fields in before and after (same number of images)
        "normalize": normalize,

        # Output
        "hdf5output": os.path.join(destpath,"{}{}.h5".format(radix[0],extout)),
        "roi": roiraw,
        "rebin":rebin}
    return config
    
def createconfig_pre(sourcepath,destpath,radix,**kwargs):
    config = createconfig(sourcepath,destpath,radix,**kwargs)

    # Create configuration file
    if not os.path.exists(destpath):
        os.makedirs(destpath)
    jsonfile = os.path.join(destpath,"{}{}.json".format(radix[0],extout))
    with open(jsonfile,'w') as f:
        json.dump(config,f,indent=2)

    return jsonfile,config["hdf5output"]

def process(sourcepath,destpath,radix,**kwargs):

    # Parse parameters
    # ... stack
    bsamefile = False
    stackdim = kwargs.get("stackdim",2)
    # ... align
    alignmethod = kwargs.get("alignmethod",None)
    refimageindex = kwargs.get("refimageindex",None)
    roialign = kwargs.get("roialign",None)
    roiresult = kwargs.get("roiresult",None)
    cropalign = kwargs.get("crop",None)
    plot = kwargs.get("plot",False)
    # ... normalization
    skipnormalization = kwargs.get("skipnormalization",False)
    normalizeonload = kwargs.get("normalizeonload",True)
    flatbefore = kwargs.get("flatbefore",True)
    flatafter = kwargs.get("flatafter",True)
    # ... other
    skippre = kwargs.get("skippre",False)
    
    logger = logging.getLogger(__name__)
    T0 = timing.taketimestamp()

    # Image stack
    normalize = normalizeonload and not skipnormalization
    jsonfile, h5file = createconfig_pre(sourcepath,destpath,radix,**kwargs)
    preprocessingexists = False
    if skippre:
        preprocessingexists = os.path.isfile(h5file)

    if preprocessingexists:
        stacks, axes, procinfo = getstacks(h5file,["detector0"])
    else:
        logger.info("Creating image stacks ...")
        stacks, axes = makestacks(jsonfile)

        #stacks2, axes2 = getstacks(h5file,["detector0"])
        #assert(axes == axes2)
        #assert(stacks == stacks2)

    h5filelast = h5file
    
    # Convert stack dictionary to stack list
    stacks = flattenstacks(stacks)

    # Default group
    default = "sample"
    defaultstack(h5file,stacks,default)

    # Groups that don't change and need to be copied
    copygroups = None

    # I0 normalization
    if not skipnormalization and not normalizeonload:
        if flatbefore and flatafter:
            if any("flat2" in s for s in stacks):
                expression = "-ln(2*{}/({flat1}+{flat2}))"
            else:
                expression = "-ln({}/{flat1})"
        elif flatbefore:
            expression = "-ln({}/{flat1})"
        elif flatafter:
            expression = "-ln({}/{flat2})"
        else:
            logger.error("Nothing to normalize with.")
            raise ValueError("Set flatbefore or flatafter to True.")
        h5file,stacks,axes = normalizefunc(h5file,stacks,axes,copygroups,bsamefile,default,expression,["flat1","flat2"],stackdim=stackdim,copyskipped=False,extension="norm")
        h5filelast = h5file
        
    # Alignment
    if alignmethod is not None and alignmethod:
        alignreference = "sample"
        h5file,stacks,axes = align(h5file,stacks,axes, copygroups, bsamefile, default,\
                                alignmethod, alignreference, refimageindex, cropalign,\
                                roialign, plot, stackdim)
        h5filelast = h5file
        
    # Crop
    if roiresult is not None:
        cropinfo = {"roi":roiresult,"stackdim":stackdim,"reference set":"sample"}
        h5file,stacks,axes = execcrop(h5file,stacks,axes, copygroups, bsamefile, default, cropinfo)
        h5filelast = h5file
        
    timing.printtimeelapsed(T0,logger)

    return h5filelast


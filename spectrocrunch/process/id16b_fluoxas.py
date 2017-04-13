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

import os
import json
import logging
import numpy as np

from spectrocrunch.xrf.create_hdf5_imagestacks import create_hdf5_imagestacks as makestacks
from spectrocrunch.h5stacks.get_hdf5_imagestacks import get_hdf5_imagestacks as getstacks
import spectrocrunch.common.timing as timing
import spectrocrunch.io.nexus as nexus

from .proc_math import execute as normalize
from .proc_align import execute as align
from .proc_replacevalue import execute as replacevalue
from .proc_crop import execute as execcrop
from . proc_common import defaultstack
from . proc_common import flattenstacks

def createconfig_pre(sourcepath,destpath,scanname,scannumbers,cfgfiles,dtcor,stackdim,\
                    counters=[],normcounter=None,mlines={},addbeforefit=True,\
                    exclude_detectors=[]):
    bfit = cfgfiles is not None

    if not isinstance(sourcepath,list):
        sourcepath = [sourcepath]
    if not isinstance(scanname,list):
        scanname = [scanname]
    if not isinstance(scannumbers,list):
        scannumbers = [scannumbers]
    if not isinstance(cfgfiles,list):
        cfgfiles = [cfgfiles]

    config = {
            # Input
            "sourcepath": sourcepath,
            "scanname": scanname,
            "scannumbers": scannumbers,
            "counters": counters,
            "counter_reldir": ".",

            # Meta data
            "metacounters": ["xia"], # get meta data from the spectrum headers
            "stacklabel": "energy",
            "scanlabel": "title",
            "coordinates": ["sy", "sz", "sx", "sampy", "sampz", "sypz", "szpz"],

            # Deadtime correction
            "dtcor": dtcor,
            "dtcorifsingle":True,

            # Configuration for fitting
            "detectorcfg": cfgfiles,
            "mlines": mlines,
            "fit": bfit,
            "fastfitting": True,
            "addbeforefitting": addbeforefit, # sum spectra
            "addafterfitting": True, # sum fit results and detector counters
            "exclude_detectors":exclude_detectors,

            # Output directories
            "outdatapath": os.path.join(destpath,scanname[0]+"_data"),
            "outfitpath": os.path.join(destpath,scanname[0]+"_fit"),
            "hdf5output": os.path.join(destpath,scanname[0]+".h5"),
            "stackdim": stackdim
    }

    # Create configuration file
    if not os.path.exists(destpath):
        os.makedirs(destpath)
    jsonfile = os.path.join(destpath,scanname[0]+".json")
    with open(jsonfile,'w') as f:
        json.dump(config,f,indent=2)

    return jsonfile,config["hdf5output"]

def process(sourcepath,destpath,scanname,scannumbers,cfgfiles,alignmethod,alignreference,\
        refimageindex=None,skippre=False,skipnormalization=False,dtcor=True,default=None,\
        crop=False,roialign=None,plot=True,counters=[],normcounter=None,mlines={},\
        addbeforefit=True,exclude_detectors=[]):

    logger = logging.getLogger(__name__)
    T0 = timing.taketimestamp()

    stackdim = 2
    bsamefile = False

    cropalign = False
    cropafter = crop
    replacenan = True

    # Image stacks
    preprocessingexists = False
    if skippre:
        h5file = os.path.join(destpath,scanname[0]+".h5")
        preprocessingexists = os.path.isfile(h5file)

    if preprocessingexists:
        stacks, axes = getstacks(h5file,["counters","^detector([0-9]+|sum)$"])
    else:
        logger.info("Creating image stacks ...")
        jsonfile, h5file = createconfig_pre(sourcepath,destpath,scanname,scannumbers,\
                                            cfgfiles,dtcor,stackdim,counters=counters,\
                                            normcounter=normcounter,mlines=mlines,\
                                            addbeforefit=addbeforefit,exclude_detectors=exclude_detectors)
        stacks, axes = makestacks(jsonfile)

        #stacks2, axes2 = getstacks(h5file,["counters","detectorsum"])
        #assert(axes == axes2)
        #assert(stacks == stacks2)

    # Convert stack dictionary to stack list
    stacks = flattenstacks(stacks)

    # Default group
    defaultstack(h5file,stacks,default)

    # Groups that don't change and need to be copied
    copygroups = ["coordinates"]

    # I0 normalization
    if not skipnormalization and normcounter is not None:
        expression = "{{}}/{{{}}}".format(normcounter)
        skip = [normcounter]
        h5file,stacks,axes = normalize(h5file,stacks,axes,copygroups,bsamefile,default,expression,skip,stackdim=stackdim,extension="norm")

    # Alignment
    if alignmethod is None or alignreference is None:
        timing.printtimeelapsed(T0,logger)
        return
    else:
        h5file,stacks,axes = align(h5file,stacks,axes, copygroups, bsamefile, default,\
            alignmethod, alignreference, refimageindex, cropalign, roialign, plot, stackdim)

    # Remove NaN's
    if replacenan:
        orgvalue = np.nan
        newvalue = 0
        replacevalue(h5file, stacks, axes, copygroups, bsamefile, default, orgvalue, newvalue)

    # Crop
    if cropafter:
        cropinfo = {"nanval":np.nan,"stackdim":stackdim,"reference set":alignreference}
        execcrop(h5file, stacks, axes, copygroups, bsamefile, default, cropinfo)

    timing.printtimeelapsed(T0,logger)


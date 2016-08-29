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

from .proc_normalize import execute as normalize
from .proc_align import execute as align
from .proc_replacevalue import execute as replacevalue
from .proc_crop import execute as execcrop
from . proc_common import defaultstack
from . proc_common import flattenstacks

def createconfig_pre(sourcepath,destpath,scanname,scannumbers,cfgfiles,dtcor,stackdim,counters=[],mlines={}):
    bfit = cfgfiles is not None

    if not isinstance(sourcepath,list):
        sourcepath = [sourcepath]
    if not isinstance(scanname,list):
        scanname = [scanname]
    if not isinstance(scannumbers,list):
        scannumbers = [scannumbers]

    if len(counters)==0:
        counters = ["zap_p201_I0","zap_p201_IC","zap_p201_It","arr_srcurr","xmap_nika","xmap_tika","xmap_feka","xmap_kka","xmap_caka","xmap_crka"]
    else:
        if "arr_srcurr" not in counters:
            counters.append("arr_srcurr")

    config = {
            # Input
            "sourcepath": sourcepath,
            "scanname": scanname,
            "scannumbers": scannumbers,
            "counters": counters,

            # Meta data
            "metacounters": ["xia"], # get meta data from the spectrum headers
            "stacklabel": "energy",
            "scanlabel": "title",
            "coordinates": ["sy", "sz", "sx", "sampy", "sampz", "sypz", "szpz"],

            # Deadtime correction
            "dtcor": dtcor,

            # Configuration for fitting
            "detectorcfg": cfgfiles,
            "mlines": mlines,
            "fit": bfit,
            "fastfitting": True,
            "addbeforefitting": False, # sum spectra
            "addafterfitting": True, # sum fit results and detector counters
            "exclude_detectors":[0,1,5],

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

def process(sourcepath,destpath,scanname,scannumbers,cfgfile,alignmethod,alignreference,\
        refimageindex=None,skippre=False,skipnormalization=False,dtcor=True,default=None,\
        crop=False,roialign=None,plot=True,counters=[],mlines={}):

    logger = logging.getLogger(__name__)
    T0 = timing.taketimestamp()

    stackdim = 2
    bsamefile = False

    cropalign = False
    cropafter = crop
    replacenan = True

    # Image stacks
    if skippre:
        h5file = os.path.join(destpath,scanname[0]+".h5")
        stacks, axes = getstacks(h5file,["counters","detectorsum"])
    else:
        logger.info("Creating image stacks ...")
        jsonfile, h5file = createconfig_pre(sourcepath,destpath,scanname,scannumbers,cfgfile,dtcor,stackdim,counters=counters,mlines=mlines)
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
    if not skipnormalization:
        h5file,stacks,axes = normalize(h5file,stacks,axes,copygroups,bsamefile,default,["arr_srcurr"],["arr_srcurr"])

    # Alignment
    if alignmethod is None or alignreference is None:
        timing.printtimeelapsed(T0,logger)
        exit()
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


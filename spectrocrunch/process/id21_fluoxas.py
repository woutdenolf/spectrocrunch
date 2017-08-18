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

from ..xrf.create_hdf5_imagestacks import create_hdf5_imagestacks as makestacks
from ..h5stacks.get_hdf5_imagestacks import get_hdf5_imagestacks as getstacks
import ..common.timing as timing
import ..io.nexus as nexus

from .proc_math import execute as normalize
from .proc_align import execute as align
from .proc_replacevalue import execute as replacevalue
from .proc_crop import execute as execcrop
from .proc_common import defaultstack
from .proc_common import flattenstacks

def createconfig_pre(sourcepath,destpath,scanname,scannumbers,cfgfile,dtcor,stackdim,mlines={},microdiff=False):
    bfit = cfgfile is not None
    if not bfit:
        cfgfile = os.path.join(destpath,scanname[0]+".cfg")

    if not isinstance(sourcepath,list):
        sourcepath = [sourcepath]
    if not isinstance(scanname,list):
        scanname = [scanname]
    if not isinstance(scannumbers,list):
        scannumbers = [scannumbers]

    if microdiff:
        counters = ["zap_iodet","zap_idet","xmap_x1","xmap_x1c","xmap_x2","xmap_x2c","xmap_x3","xmap_x3c","xmap_icr","xmap_ocr"]
        motors = ["samh", "samv", "samd", "samph", "sampv"]
        counter_reldir = ".."
    else:
        counters = ["arr_iodet","arr_idet","arr_absorp1","arr_absorp2","arr_absorp3","xmap_x1","xmap_x1c","xmap_x2","xmap_x2c","xmap_x3","xmap_x3c","xmap_icr","xmap_ocr"]
        motors = ["samy", "samz", "samx", "sampy", "sampz"]
        counter_reldir = "."

    config = {
            # Input
            "sourcepath": sourcepath,
            "scanname": scanname,
            "scannumbers": scannumbers,
            "counters": counters,
            "counter_reldir": counter_reldir,

            # Meta data
            "metacounters": counters,
            "stacklabel": "DCM_Energy",
            "fastlabel": "fast",
            "slowlabel": "slow",
            "coordinates": motors,

            # Deadtime correction
            "dtcor": dtcor,

            # Configuration for fitting
            "detectorcfg": [cfgfile],
            "mlines": mlines,
            "fit": bfit,
            "fastfitting": True,
            "addbeforefitting": False,
            "addafterfitting": False,
            "exclude_detectors":[],

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
        crop=False,roialign=None,plot=True,mlines={},microdiff=False):

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
        stacks, axes = getstacks(h5file,["counters","detector0"])
    else:
        logger.info("Creating image stacks ...")
        jsonfile, h5file = createconfig_pre(sourcepath,destpath,scanname,scannumbers,cfgfile,False,stackdim,mlines=mlines,microdiff=microdiff)
        stacks, axes = makestacks(jsonfile)

        #stacks2, axes2 = getstacks(h5file,["counters","detector0"])
        #assert(axes == axes2)
        #assert(stacks == stacks2)

    # Convert stack dictionary to stack list
    stacks = flattenstacks(stacks)

    # Default group
    defaultstack(h5file,stacks,default)

    # Groups that don't change and need to be copied
    copygroups = ["coordinates"]

    # I0 normalization
    if not skipnormalization or dtcor:
        skip = ["arr_absorp1","arr_absorp2","arr_absorp3"]
        if microdiff:
            normcounter = "zap_iodet"
        else:
            normcounter = "arr_iodet"

        if dtcor:
            expression = "{{}}*{{xmap_icr}}/({{{}}}*{{xmap_ocr}})".format(normcounter)
            skip += [normcounter,"xmap_icr","xmap_ocr"]
        else:
            expression = "{{}}/{{{}}}".format(normcounter)
            skip += [normcounter]

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

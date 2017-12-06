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
import copy

#from ..xrf.create_hdf5_imagestacks_old import create_hdf5_imagestacks as makestacks
from ..xrf.create_hdf5_imagestacks import create_hdf5_imagestacks as makestacks
from ..h5stacks.get_hdf5_imagestacks import get_hdf5_imagestacks as getstacks
from ..common import timing
from ..io import nexus

from .proc_math import execute as math
from .proc_align import execute as align
from .proc_replacevalue import execute as replacevalue
from .proc_crop import execute as execcrop
from .proc_common import defaultstack
from .proc_common import flattenstacks
from .proc_resample import execute as execresample

logger = logging.getLogger(__name__)

def createconfig_pre(sourcepath,destpath,scanname,scannumbers,cfgfiles,dtcor,stackdim,\
                    mlines={},microdiff=False,exclude_detectors=None,include_detectors=None,\
                    addbeforefit=True,useencoders=False,noxia=False):
    if noxia:
        cfgfiles = None
    bfit = cfgfiles is not None

    if exclude_detectors is None:
        exclude_detectors = []

    if include_detectors is None:
        include_detectors = []

    if not isinstance(sourcepath,list):
        sourcepath = [sourcepath]
    if not isinstance(scanname,list):
        scanname = [scanname]
    if not isinstance(scannumbers,list):
        scannumbers = [scannumbers]
    if not isinstance(cfgfiles,list):
        cfgfiles = [cfgfiles]

    if microdiff:
        if noxia:
            counters = ["zap_iodet","zap_idet"]
        else:
            counters = ["zap_iodet","zap_idet","xmap_x1","xmap_x1c","xmap_x2","xmap_x2c","xmap_x3","xmap_x3c","xmap_icr","xmap_ocr"]
        motors = ["samh", "samv", "samd", "samph", "sampv"]
        counter_reldir = ".."
    else:
        if noxia:
            counters = ["arr_iodet","arr_idet","arr_fdet","arr_absorp1","arr_absorp2","arr_absorp3"]
        else:
            counters = ["arr_iodet","arr_idet","arr_fdet","arr_absorp1","arr_absorp2","arr_absorp3","xmap_x1","xmap_x1c","xmap_x2","xmap_x2c","xmap_x3","xmap_x3c","xmap_icr","xmap_ocr"]
        if useencoders:
            counters += ["arr_samy","arr_samz"]
        motors = ["samy", "samz", "samx", "sampy", "sampz"]
        counter_reldir = "."

    config = {
            # Input
            "sourcepath": sourcepath,
            "scanname": scanname,
            "scannumbers": scannumbers,
            "counters": counters,
            "fluxcounter": 0,
            "counter_reldir": counter_reldir,

            # Meta data
            "metacounters": counters,
            "stacklabel": "DCM_Energy",
            "speccmdlabel": "title",
            #"fastlabel": "fast",
            #"slowlabel": "slow",
            #"timelabel": "time",
            "coordinates": motors,

            # Deadtime correction
            "dtcor": dtcor,
            "dtcorifsingle":False, # correct for deadtime when a single detector (ignored when dtcor==False)

            # Configuration for fitting
            "detectorcfg": cfgfiles,
            "mlines": mlines,
            "fit": bfit,
            "fastfitting": True,
            "addbeforefitting": addbeforefit, # sum spectra
            "addafterfitting": True, # sum fit results and detector counters
            "exclude_detectors": exclude_detectors,
            "include_detectors": include_detectors,
            
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

def process(sourcepath,destpath,scanname,scannumbers,cfgfiles,\
        skippre=False,dtcor=True,exclude_detectors=None,include_detectors=None,addbeforefit=True,encodercor=None,noxia=False,mlines={},\
        default=None,microdiff=False,fluxmonitor=None,\
        alignmethod=None,alignreference=None,refimageindex=None,roialign=None,plot=True,\
        crop=False,replacenan=True,prealignnormcounter=None,postalignnormcounter=None,stackdim=2):

    T0 = timing.taketimestamp()

    bsamefile = False
    cropalign = False
    cropafter = crop
    
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
                                            cfgfiles,dtcor,stackdim,mlines=mlines,microdiff=microdiff,\
                                            exclude_detectors=exclude_detectors,include_detectors=include_detectors,addbeforefit=addbeforefit,\
                                            useencoders=encodercor is not None,noxia=noxia)
        stacks, axes = makestacks(jsonfile,fluxmonitor=fluxmonitor)

        #stacks2, axes2 = getstacks(h5file,["counters","detectorsum"])
        #assert(axes == axes2)
        #assert(stacks == stacks2)
        
        
    if "detectorsum" in stacks:
        dtcor = False # done on the raw data

    # Convert stack dictionary to stack list
    stacks = flattenstacks(stacks)

    # Default group
    defaultstack(h5file,stacks,default)

    # Groups that don't change and need to be copied
    copygroups = ["coordinates"]

    # Normalization
    skipnorm = ["arr_absorp1","arr_absorp2","arr_absorp3","arr_samy","arr_samz"]
    if dtcor or prealignnormcounter is not None:
        skip = copy.copy(skipnorm)

        # Create normalization expression
        if dtcor:
            if prealignnormcounter is None:
                expression = "{{}}*nanone({{xmap_icr}}/{{xmap_ocr}})"
            else:
                expression = "{{}}*nanone({{xmap_icr}}/({{{}}}*{{xmap_ocr}}))".format(prealignnormcounter)
            skip += ["xmap_icr","xmap_ocr"]
        else:
            expression = "{{}}/{{{}}}".format(prealignnormcounter)

        if prealignnormcounter is not None:
            skip += [prealignnormcounter]

        h5file,stacks,axes = math(h5file,stacks,axes,copygroups,bsamefile,default,expression,skip,stackdim=stackdim,extension="norm")

    # Correct for encoder positions
    if encodercor is not None:
        resampleinfo = {}

        if "samy" in encodercor:
            resampleinfo["samy"] = {"encoder":"arr_samy","resolution":encodercor["samy"]} # resolution: 52500 steps/unit
            
        if "samz" in encodercor:
            resampleinfo["samz"] = {"encoder":"arr_samz","resolution":encodercor["samz"]} # resolution: 50000 steps/unit

        if len(resampleinfo) is not None:
            h5file,stacks,axes = execresample(h5file, stacks, axes, copygroups, bsamefile, default, resampleinfo)

    # Alignment
    if alignmethod is not None and alignreference is not None:
        # Alignment
        h5file,stacks,axes = align(h5file,stacks,axes, copygroups, bsamefile, default,\
            alignmethod, alignreference, refimageindex, cropalign, roialign, plot, stackdim)

        # Post normalization
        if postalignnormcounter is not None:
            skip = skipnorm+[postalignnormcounter]
            expression = "{{}}/{{{}}}".format(postalignnormcounter)
            h5file,stacks,axes = math(h5file,stacks,axes,copygroups,bsamefile,default,expression,skip,stackdim=stackdim,extension="postnorm")

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

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

from spectrocrunch.xrf.create_hdf5_imagestacks import create_hdf5_imagestacks as makestacks
from spectrocrunch.xrf.get_hdf5_imagestacks import get_hdf5_imagestacks as getstacks
from spectrocrunch.xrf.math_hdf5_imagestacks import fluxnorm_hdf5_imagestacks as fluxnormstacks
from spectrocrunch.xrf.math_hdf5_imagestacks import copy_hdf5_imagestacks as copystacks
from spectrocrunch.xrf.align_hdf5_imagestacks import align_hdf5_imagestacks as alignstacks
import spectrocrunch.common.timing as timing
import spectrocrunch.io.nexus as nexus

def createconfig_pre(sourcepath,destpath,scanname,scannumbers,cfgfile,dtcor,stackdim):
    bfit = cfgfile is not None
    if not bfit:
        cfgfile = os.path.join(destpath,scanname[0]+".cfg")

    if not isinstance(sourcepath,list):
        sourcepath = [sourcepath]
    if not isinstance(scanname,list):
        scanname = [scanname]
    if not isinstance(scannumbers,list):
        scannumbers = [scannumbers]

    config = {
            # Input
            "sourcepath": sourcepath,
            "scanname": scanname,
            "scannumbers": scannumbers,
            "counters": ["arr_iodet","arr_idet","arr_absorp1","arr_absorp2","arr_absorp3","xmap_x1","xmap_x1c","xmap_x2","xmap_x2c","xmap_x3","xmap_x3c"],

            # Meta data
            "metacounters": ["arr_iodet","arr_idet","arr_absorp1","arr_absorp2","arr_absorp3","xmap_x1","xmap_x1c","xmap_x2","xmap_x2c","xmap_x3","xmap_x3c"],
            "stacklabel": "DCM_Energy",
            "fastlabel": "fast",
            "slowlabel": "slow",
            "coordinates": ["samy", "samz", "samx", "sampy", "sampz"],

            # Deadtime correction
            "dtcor": dtcor,

            # Configuration for fitting
            "detectorcfg": [cfgfile],
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

def isgroup(name,reference):
    return name.split(".")[0].endswith(reference)

def defaultstack(f,stacks,plotreference):
    if plotreference is None:
        reference = []
    else:
        reference = [s for s in stacks if isgroup(s,plotreference)]
    if len(reference)!=1:
        return
    nexus.defaultstack(f,reference[0])

def process(sourcepath,destpath,scanname,scannumbers,cfgfile,alignmethod,alignreference,refimageindex=None,skippre=False,skipnormalization=False,dtcor=True,default=None,crop=False):

    logger = logging.getLogger(__name__)
    T0 = timing.taketimestamp()

    stackdim = 2
    bsamefile = False

    # Image stack
    logger.info("Creating image stacks ...")
    if skippre:
        filein = os.path.join(destpath,scanname[0]+".h5")
        stacks, axes = getstacks(filein,["counters","detector0"])
    else:
        jsonfile, filein = createconfig_pre(sourcepath,destpath,scanname,scannumbers,cfgfile,dtcor,stackdim)
        stacks, axes = makestacks(jsonfile)

        #stacks2, axes2 = getstacks(filein,["counters","detector0"])
        #assert(axes == axes2)
        #assert(stacks == stacks2)

    # Default
    defaultstack(filein,stacks["counters"].values(),default)
    if "detector0" in stacks:
        defaultstack(filein,stacks["detector0"].values(),default)

    # Groups that don't change and need to be copied
    copygroups = ["coordinates"]

    # I0 normalization
    base,ext = os.path.splitext(filein)
    if "detector0" in stacks and "arr_iodet" in stacks["counters"] and skipnormalization==False:
        logger.info("I0 normalization ...")
        if bsamefile:
            fileout = filein
        else:
            fileout = base+".norm"+ext

        # normalization stacks
        I0stack = stacks["counters"]["arr_iodet"]

        # stacks to be normalized
        innames = stacks["detector0"].values()
        notnormalize = ["arr_absorp1","arr_absorp2","arr_absorp3"]
        for ctr in stacks["counters"]:
            if not any(isgroup(ctr,s) for s in notnormalize):
                innames += [stacks["counters"][ctr]]

        # normalize
        Ifn_stacks, Ifn_axes = fluxnormstacks(filein,fileout,axes,I0stack,innames,innames,overwrite=True,info={"normalization":"arr_iodet"},copygroups=copygroups)

        # copy unnormalized stacks when new file
        innames = [ctr for ctr in stacks["counters"].values() if any(isgroup(ctr,s) for s in notnormalize)]
        if filein==fileout:
            Ifn_stacks += innames
        else:
            tmp_stacks, tmp = copystacks(filein,fileout,axes,innames,innames,overwrite=False)
            Ifn_stacks += tmp_stacks

        filein = fileout

        # Default
        defaultstack(fileout,Ifn_stacks,default)

    else:
        Ifn_stacks = stacks["counters"].values()
        if "detector0" in stacks:
            Ifn_stacks += [s for s in stacks["detector0"].values()]
        Ifn_axes = [dict(a) for a in axes]

    # Alignment
    if alignmethod is not None and alignreference is not None:
        logger.info("Aligning image stacks ...")
        if bsamefile:
            fileout = filein
        else:
            fileout = base+".align"+ext

        info = {"method":alignmethod,"pairwise":refimageindex==None,\
                "reference set":alignreference,\
                "reference image":"None" if refimageindex==None else refimageindex,\
                "crop":crop}
        aligned_stacks, aligned_axes = alignstacks(filein,Ifn_stacks,Ifn_axes,stackdim,fileout,alignmethod,
                                        alignreference,refimageindex=refimageindex,overwrite=True,crop=crop,
                                        info=info,copygroups=copygroups)
        # Default
        defaultstack(fileout,aligned_stacks,default)

    timing.printtimeelapsed(T0,logger)


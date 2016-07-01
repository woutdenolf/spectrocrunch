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
from spectrocrunch.h5stacks.math_hdf5_imagestacks import fluxnorm_hdf5_imagestacks as fluxnormstacks
from spectrocrunch.h5stacks.math_hdf5_imagestacks import copy_hdf5_imagestacks as copystacks
from spectrocrunch.h5stacks.math_hdf5_imagestacks import crop_hdf5_imagestacks as cropstacks
from spectrocrunch.h5stacks.math_hdf5_imagestacks import replacevalue_hdf5_imagestacks as replacestacks
from spectrocrunch.h5stacks.align_hdf5_imagestacks import align_hdf5_imagestacks as alignstacks
import spectrocrunch.common.timing as timing
import spectrocrunch.io.nexus as nexus

def createconfig_pre(sourcepath,destpath,scanname,scannumbers,cfgfiles,dtcor,stackdim):
    bfit = cfgfiles is not None

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
            "counters": ["zap_norm_It","zap_p201_I0","zap_p201_IC","zap_p201_It","arr_srcurr","xmap_tika","xmap_ska","xmap_feka","xmap_auka","xmap_aska"],

            # Meta data
            "metacounters": ["xia"], # get meta data from the spectrum headers
            "stacklabel": "energy",
            "scanlabel": "title",
            "coordinates": ["sy", "sz", "sx", "sampy", "sampz", "sypz", "szpz"],

            # Deadtime correction
            "dtcor": dtcor,

            # Configuration for fitting
            "detectorcfg": cfgfiles,
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

def process(sourcepath,destpath,scanname,scannumbers,cfgfile,alignmethod,alignreference,\
        refimageindex=None,skippre=False,skipnormalization=False,dtcor=True,default=None,crop=False,roi=None,plot=True):

    logger = logging.getLogger(__name__)
    T0 = timing.taketimestamp()

    stackdim = 2
    bsamefile = False

    cropalign = crop
    cropafter = False
    replacenan = False

    # Image stack
    logger.info("Creating image stacks ...")
    if skippre:
        file_raw = os.path.join(destpath,scanname[0]+".h5")
        stacks, axes = getstacks(file_raw,["counters","detectorsum"])
    else:
        jsonfile, file_raw = createconfig_pre(sourcepath,destpath,scanname,scannumbers,cfgfile,dtcor,stackdim)
        stacks, axes = makestacks(jsonfile)

        #stacks2, axes2 = getstacks(file_raw,["counters","detectorsum"])
        #assert(axes == axes2)
        #assert(stacks == stacks2)

    # Default
    defaultstack(file_raw,stacks["counters"].values(),default)
    if "detectorsum" in stacks:
        defaultstack(file_raw,stacks["detectorsum"].values(),default)

    # Groups that don't change and need to be copied
    copygroups = ["coordinates"]

    # I0 normalization and convert stack dictionary to stack list
    base,ext = os.path.splitext(file_raw)
    if "detectorsum" in stacks and "arr_srcurr" in stacks["counters"] and skipnormalization==False:
        logger.info("I0 normalization ...")
        if bsamefile:
            file_normalized = file_raw
        else:
            file_normalized = base+".norm"+ext

        # normalization stacks
        I0stacks = [stacks["counters"]["arr_srcurr"]]

        # stacks to be normalized
        innames = stacks["detectorsum"].values()
        notnormalize = []
        for ctr in stacks["counters"]:
            if not any(isgroup(ctr,s) for s in notnormalize):
                innames += [stacks["counters"][ctr]]

        # normalize
        Ifn_stacks, Ifn_axes = fluxnormstacks(file_raw,file_normalized,axes,I0stacks,innames,innames,overwrite=True,info={"normalization":"arr_srcurr"},copygroups=copygroups)

        # copy unnormalized stacks when new file
        innames = [ctr for ctr in stacks["counters"].values() if any(isgroup(ctr,s) for s in notnormalize)]
        if file_raw==file_normalized:
            Ifn_stacks += innames
        else:
            tmp_stacks, tmp = copystacks(file_raw,file_normalized,axes,innames,innames,overwrite=False)
            Ifn_stacks += tmp_stacks

        # Default
        defaultstack(file_normalized,Ifn_stacks,default)
    else:
        Ifn_stacks = stacks["counters"].values()
        if "detectorsum" in stacks:
            Ifn_stacks += [s for s in stacks["detectorsum"].values()]
        Ifn_axes = [dict(a) for a in axes]
        file_normalized = file_raw

    # Alignment
    if alignmethod is not None and alignreference is not None:
        logger.info("Aligning image stacks ...")
        if bsamefile:
            file_aligned = file_normalized
        else:
            file_aligned = base+".align"+ext

        info = {"method":alignmethod,"pairwise":refimageindex==None,\
                "reference set":alignreference,\
                "reference image":refimageindex,\
                "crop":cropalign,\
                "roi":roi}
        aligned_stacks, aligned_axes = alignstacks(file_normalized,Ifn_stacks,Ifn_axes,stackdim,file_aligned,alignmethod,
                                        alignreference,refimageindex=refimageindex,overwrite=True,crop=cropalign,
                                        roi=roi,plot=plot,info=info,copygroups=copygroups)

        # Default
        defaultstack(file_aligned,aligned_stacks,default)
    else:
        aligned_stacks = Ifn_stacks
        aligned_axes = Ifn_axes
        file_aligned = file_normalized

    # Remove NaN's
    if replacenan and alignmethod is not None and alignreference is not None:
        logger.info("Replace NaN's ...")
        if bsamefile:
            file_nonan = file_aligned
        else:
            file_nonan = base+".nonan"+ext

        orgvalue = np.nan
        newvalue = 0
        info = {"replaced value":orgvalue,"new value":newvalue}
        replaced_stacks, replaced_axes = replacestacks(file_aligned,file_nonan,aligned_axes,aligned_stacks,aligned_stacks,orgvalue,newvalue,overwrite=True,info=info,copygroups=copygroups)

        # Default
        defaultstack(file_nonan,replaced_stacks,default)
    else:
        replaced_stacks = aligned_stacks
        replaced_axes = aligned_axes
        file_nonan = file_aligned

    # Crop
    if cropafter and alignreference is not None:
        logger.info("Crop image stacks ...")
        if bsamefile:
            file_cropped = file_aligned
        else:
            file_cropped = base+".crop"+ext

        info = {"crop value":np.nan,"reference set":alignreference}
        cropinfo = {"nanval":np.nan,"stackdim":stackdim,"reference set":alignreference}
        cropped_stacks, cropped_axes = cropstacks(file_aligned,file_cropped,aligned_axes,aligned_stacks,aligned_stacks,cropinfo,overwrite=True,info=info,copygroups=copygroups)

        # Default
        defaultstack(file_cropped,cropped_stacks,default)
    else:
        cropped_stacks = aligned_stacks
        cropped_axes = aligned_axes
        file_cropped = file_aligned
    
    timing.printtimeelapsed(T0,logger)


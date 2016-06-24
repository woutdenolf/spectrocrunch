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
from spectrocrunch.h5stacks.math_hdf5_imagestacks import fluxnorm_hdf5_imagestacks as fluxnormstacks
from spectrocrunch.h5stacks.math_hdf5_imagestacks import copy_hdf5_imagestacks as copystacks
from spectrocrunch.h5stacks.math_hdf5_imagestacks import crop_hdf5_imagestacks as cropstacks
from spectrocrunch.h5stacks.align_hdf5_imagestacks import align_hdf5_imagestacks as alignstacks
import spectrocrunch.io.nexus as nexus

def createconfig_pre(sourcepath,destpath,radix,ext,rebin,stackdim):

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
        "dtype": "np.float32",
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
        refimageindex=None,crop=False,roialign=None,roiresult=None,plot=True):

    logger = logging.getLogger(__name__)
    T0 = timing.taketimestamp()

    stackdim = 2
    bsamefile = False

    # Image stack
    jsonfile, file_raw = createconfig_pre(sourcepath,destpath,radix,ext,rebin,stackdim)
    if skippre:
        stacks, axes = getstacks(file_raw,["detector0"])
    else:
        logger.info("Creating image stacks ...")
        stacks, axes = makestacks(jsonfile)
        nexus.defaultstack(file_raw,stacks["detector0"]["sample"])

        #stacks2, axes2 = getstacks(file_raw,["detector0"])
        #assert(axes == axes2)
        #assert(stacks == stacks2)

    # I0 normalization and convert stack dictionary to stack list
    base,ext = os.path.splitext(file_raw)
    if bsamefile:
        file_normalized = file_raw
    else:
        file_normalized = base+".norm"+ext
    if skipnormalization:
        Ifn_stacks, Ifn_axes = getstacks(file_normalized,["detector0"])
        Ifn_stacks = Ifn_stacks["detector0"].values()
    else:
        logger.info("I0 normalization ...")
        # normalization stacks
        if "flat2" in stacks["detector0"]:
            I0stacks = [stacks["detector0"]["flat1"],stacks["detector0"]["flat2"]]
            info = {"normalization":"(flat1+flat2)/2"}
        else:
            I0stacks = [stacks["detector0"]["flat1"]]
            info = {"normalization":"flat1"}

        # stacks to be normalized
        innames = [stacks["detector0"]["sample"]]

        # normalize
        Ifn_stacks, Ifn_axes = fluxnormstacks(file_raw,file_normalized,axes,I0stacks,innames,innames,overwrite=True,info=info,stackdim=stackdim)

        # copy unnormalized stacks when new file
        #innames = I0stacks
        #if file_raw==file_normalized:
        #    Ifn_stacks += innames
        #else:
        #    tmp_stacks, tmp = copystacks(file_raw,file_normalized,axes,innames,innames,overwrite=False)
        #    Ifn_stacks += tmp_stacks

        # Default
        nexus.defaultstack(file_normalized,Ifn_stacks[0])

    # Alignment
    if bsamefile:
        file_aligned = file_normalized
    else:
        file_aligned = base+".align"+ext
    if skipalign:
        aligned_stacks, aligned_axes = getstacks(file_normalized,["detector0"])
        aligned_stacks = aligned_stacks["detector0"].values()
    else:
        if alignmethod is not None:
            logger.info("Aligning image stack ...")
            alignreference = "sample"
            info = {"method":alignmethod,"pairwise":refimageindex==None,\
                    "reference set":alignreference,\
                    "reference image":refimageindex,\
                    "crop":crop,\
                    "roi":roialign}
            aligned_stacks, aligned_axes = alignstacks(file_normalized,Ifn_stacks,Ifn_axes,stackdim,file_aligned,alignmethod,
                                            alignreference,refimageindex=refimageindex,overwrite=True,crop=crop,
                                            roi=roialign,plot=plot,info=info)

            # Default
            nexus.defaultstack(file_aligned,aligned_stacks[0])
        else:
            aligned_stacks = Ifn_stacks
            aligned_axes = Ifn_axes
            file_aligned = file_normalized

    # Crop result
    if roiresult is not None:
        logger.info("Crop result ...")

        if bsamefile:
            file_cropped = file_aligned
        else:
            file_cropped = base+".crop"+ext

        info = {"image roi":roiresult}
        cropinfo = {"roi":roiresult,"stackdim":stackdim,"reference set":"sample"}
        cropped_stacks, cropped_axes = cropstacks(file_aligned,file_cropped,aligned_axes,aligned_stacks,aligned_stacks,cropinfo,overwrite=True,info=info)

        # Default
        nexus.defaultstack(file_cropped,cropped_stacks[0])

    timing.printtimeelapsed(T0,logger)



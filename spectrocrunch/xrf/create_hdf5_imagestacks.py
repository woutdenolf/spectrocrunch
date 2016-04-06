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
import h5py
from PyMca5.PyMcaIO import EdfFile

from spectrocrunch.xrf.parse_xia import parse_xia_esrf
from spectrocrunch.xrf.fit import PerformBatchFit as fitter
from . import preparenexus_hdf5_imagestacks as nexus

def filecounter(sourcepath,scanname,counter,scannumber,idet=None):
    if idet is not None:
        filename = "%s_%s_%02d_%04d_0000.edf"%(scanname,counter,idet,scannumber)
    else:
        filename = "%s_%s_%04d_0000.edf"%(scanname,counter,scannumber)
    counterfile = os.path.join(sourcepath,filename)
    return counterfile

def detectorname(config,ndet,idet,counter=False):
    if ndet == 0:
        return ""
    elif config["addbeforefitting"] and ndet > 1 and not counter:
        return "sum"
    else:
        return "detector%d"%idet

def getimagestacks(config):
    """Get image stacks (counters, ROI's, fitted maps)
    """

    # Prepare data
    npaths = len(config["sourcepath"])
    if npaths != len(config["scanname"]):
        raise ValueError("Number or scan names must be the same as number of source paths.")
    if npaths != len(config["scannumbers"]):
        raise ValueError("Number or scan number sets must be the same as number of source paths.")

    stacks = {}
    auxstacks = {}
    stackaxes = [None]*3
    coordinates = []
    nscanstot = 0
    iscanoffset = 0
    for ipath in range(npaths):
        nscanstot += len(config["scannumbers"][ipath])

    # DT correction, sum detector and fit
    nfilestofit = None
    for ipath in range(npaths):
        sourcepath = config["sourcepath"][ipath]
        scanname = config["scanname"][ipath]
        nscans = len(config["scannumbers"][ipath])

        for iscan in range(nscans):
            scannumber = config["scannumbers"][ipath][iscan]

            # Extract stack axes values
            stackvalue = np.nan
            for metacounter in config["metacounters"]:
                try:
                    metafile = filecounter(sourcepath,scanname,metacounter,scannumber,idet=0 if "xmap" in metacounter else None)
                    print(metafile)
                    metafile = EdfFile.EdfFile(metafile)
                    

                    if iscan == 0 and ipath == 0:
                        label = config["fastlabel"]
                        motfast = metafile.Images[0].Header[label+"_mot"]
                        start = np.float(metafile.Images[0].Header[label+"_start"])
                        end = np.float(metafile.Images[0].Header[label+"_end"])
                        nbp = np.float(metafile.Images[0].Header[label+"_nbp"])
                        stackaxes[1] = {"name":str(motfast),"data":np.linspace(start,end,nbp)}

                        label = config["slowlabel"]
                        motslow = metafile.Images[0].Header[label+"_mot"]
                        start = np.float(metafile.Images[0].Header[label+"_start"])
                        end = np.float(metafile.Images[0].Header[label+"_end"])
                        nbp = np.float(metafile.Images[0].Header[label+"_nbp"])
                        end += (end-start)/(nbp-1)
                        nbp += 1
                        stackaxes[0] = {"name":str(motslow),"data":np.linspace(start,end,nbp)}

                        stackaxes[2] = {"name":str(config["stacklabel"]),"data":np.full(nscanstot,np.nan,dtype=np.float32)}

                        coordinates = {mot:np.float32(metafile.Images[0].Header[mot]) for mot in config["coordinates"] if mot != motfast and mot != motslow}

                    stackvalue = np.float(metafile.Images[0].Header[config["stacklabel"]])
                    break
                except:
                    logging.exception("Something wrong with extracting info from metacounter.")
            if stackaxes[2] is None:
                raise IOError("Metacounter files are not present, corrupted or not the right format.")
            stackaxes[2]["data"][iscanoffset+iscan] = stackvalue

            # DT correction and sum detector (if more than 1)
            if config["dtcor"]:
                parsename = "dtcor"
            else:
                parsename = "copy"
            parsename = "%%0%dd_%s"%(np.int(np.floor(np.log10(npaths)))+1,parsename)%(ipath)

            filestofit,detnums = parse_xia_esrf(sourcepath,scanname,
                    scannumber,config["outdatapath"],parsename,add=config["addbeforefitting"],
                    exclude_detectors=config["exclude_detectors"],deadtime=config["dtcor"])
            ndet = len(detnums)

            # Allocate first level of the stacks
            if iscan == 0 and ipath == 0:
                counters = [ctr for ctr in config["counters"] if "xmap" not in ctr]
                detcounters = [ctr for ctr in config["counters"] if "xmap" in ctr]

                for idet in detnums:
                    stacks[detectorname(config,ndet,idet)] = {}
                    if len(detcounters)>0:
                        stacks[detectorname(config,ndet,idet,counter=True)] = {}

                stacks["counters"] = {}
                for counter in counters:
                    stacks["counters"][counter] = [""]*nscanstot

            # Loop over the detectors
            for i in range(ndet):
                idet = detnums[i]

                if len(filestofit[idet])!= 0 and config["fit"]:
                    if len(config["detectorcfg"])==1:
                        cfg = config["detectorcfg"][0]
                    else:
                        cfg = config["detectorcfg"][i]

                    # Fit spectra resulting in images
                    outname = "%s_%s_xia%02d_%04d_0000"%(scanname,parsename,idet,scannumber)
                    files, labels = fitter(filestofit[idet],
                                                config["outfitpath"],outname,cfg,stackvalue,
                                                fast=config["fastfitting"])

                    # Append images
                    detname = detectorname(config,ndet,idet)
                    if iscan == 0 and ipath == 0:
                        for label in labels:
                            stacks[detname][label] = [""]*nscanstot
                    for i in range(len(labels)):
                        stacks[detname][labels[i]][iscanoffset+iscan] = files[i]

                # Append counters
                detname = detectorname(config,ndet,idet,counter=True)
                if iscan == 0 and ipath == 0:
                    for counter in detcounters:
                        stacks[detname][counter] = [""]*nscanstot
                for counter in detcounters:
                    stacks[detname][counter][iscanoffset+iscan] = filecounter(sourcepath,scanname,counter,scannumber,idet=idet)

            # Append counters
            for counter in counters:
                stacks["counters"][counter][iscanoffset+iscan] = filecounter(sourcepath,scanname,counter,scannumber)
        
        iscanoffset += nscans

    return stacks,stackaxes,coordinates

def exportgroups(f,stacks,keys,axes,sumgroups=False):
    """Export groups of EDF stacks, summated or not
    """

    if sumgroups:
        if "sum" in f:
            grp = f["sum"]
        else:
            grp = f.create_group("sum")

    for k1 in keys: # detector or counter group
        if not sumgroups:
            grp = f.create_group(k1)

        for k2 in stacks[k1].keys(): # stack subgroup (NXdata class)
            nxdatagrp = None

            nscans = len(stacks[k1][k2])
            for iscan in range(nscans):
                fdata = EdfFile.EdfFile(stacks[k1][k2][iscan])
                if fdata.GetNumImages() <= 0:
                    continue

                data = fdata.GetData(0)
                if k2 not in grp:
                    nxdatagrp = grp.create_group(k2)
                    nxdatagrp.attrs["NX_class"] = "NXdata"
                    nxdatagrp.attrs["signal"] = "data"
                    nxdatagrp.attrs["axes"] = [a["name"] for a in axes]

                    shape = data.shape
                    dset = nxdatagrp.create_dataset("data",shape+(nscans,),chunks = True,dtype = np.float32)
                    for a in axes:
                        nxdatagrp[a["name"]] = f[a["fullname"]]
                    if sumgroups:
                        dset[:] = 0
                    else:
                        dset[:] = np.nan

                # Some rows too much or rows missing:
                if shape[0] > data.shape[0]:
                    data = np.pad(data,((0,shape[0]-data.shape[0]),(0,0)),'constant',constant_values=self.cval)
                elif shape[0] < data.shape[0]:
                    data = data[0:shape[0],:]

                if sumgroups:
                    dset[...,iscan] += data
                else:
                    dset[...,iscan] = data

            # Replace the list-of-files by the NX entry
            if nxdatagrp is None:
                del stacks[k1][k2]
            else:
                stacks[k1][k2] = nxdatagrp.name

def exportimagestacks(config,stacks,stackaxes,coordinates,jsonfile):
    """Export EDF stacks to HDF5
    """
    f = h5py.File(config["hdf5output"],'w')

    # Save stack axes values
    grp = f.create_group("axes")
    naxes = len(stackaxes)
    axes = [None]*naxes
    for i in range(naxes):
        dset = grp.create_dataset(stackaxes[i]["name"],data=stackaxes[i]["data"])
        axes[i] = {"name":stackaxes[i]["name"],"fullname":dset.name} 

    # Save sum of detectors if requested
    keys = [k for k in stacks.keys() if "detector" in k]
    if config["addafterfitting"] and len(keys)>1:
        exportgroups(f,stacks,keys,axes,sumgroups=True)
        kleft = [k for k in stacks.keys() if "detector" not in k]
    else:
        kleft = stacks.keys()

    # Save other groups (not included in the sum)
    if len(kleft) > 0:
        exportgroups(f,stacks,kleft,axes)

    # Save coordinates
    coordgrp = f.create_group("coordinates")
    for k in coordinates:
        coordgrp[k] = coordinates[k]

    # Add processing info
    #nexus.addinfogroup(f,"fromraw",config)
    nexus.addinfogroup(f,"fromraw",{"config":jsonfile})

    f.close()

    return axes

def create_hdf5_imagestacks(jsonfile):
    """Convert scanning data (XIA spectra + counters) to an HDF5 file:
        groups which contain NXdata classes
        3 axes datasets on the main level
    """
    # Processing configuration
    with open(jsonfile,'r') as f:
        config = json.load(f)

    # Raw or pre-processed data (e.g. DT correction, fitting)
    stacks,stackaxes,coordinates = getimagestacks(config)

    # Export EDF stacks to HDF5
    axes = exportimagestacks(config,stacks,stackaxes,coordinates,jsonfile)

    return stacks,axes

if __name__ == '__main__':
    import sys
    if len(sys.argv)>=2:
        process_esrf(sys.argv[1])



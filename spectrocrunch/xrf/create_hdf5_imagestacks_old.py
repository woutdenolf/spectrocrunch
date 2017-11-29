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
from glob import glob
import fabio

from ..xrf.parse_xia import parse_xia_esrf
from ..xrf.fit import PerformBatchFit as fitter
from ..io import nexus
from ..io import spec

logger = logging.getLogger(__name__)

def filecounter_nocheck(sourcepath,scanname,counter,scannumber,idet=None):
    f = '0000'
    if counter=="xia":
        filename = "%s_%sst_%04d_0000_%s.edf"%(scanname,counter,scannumber,f)
    elif idet is not None:
        filename = "%s_%s_%02d_%04d_%s.edf"%(scanname,counter,idet,scannumber,f)
    else:
        filename = "%s_%s_%04d_%s.edf"%(scanname,counter,scannumber,f)
    return os.path.join(sourcepath,filename)
          
def filecounter(sourcepath,scanname,counter,scannumber,idet=None,getcount=False):
    if getcount:
        f = '*'
    else:
        f = '0000'

    if counter=="xia":
        filename = ["%s_%sst_%04d_0000_%s.edf"%(scanname,counter,scannumber,f),\
                    "%s_%sst_%04d_0000_%s.edf.1"%(scanname,counter,scannumber,f)]
    elif idet is not None:
        filename = ["%s_%s_%02d_%04d_%s.edf"%(scanname,counter,idet,scannumber,f),\
                    "%s_PUZ_%s_%02d_%04d_%s.edf"%(scanname,counter,idet,scannumber,f),\
                    "%s_%s_%02d_%04d_%s.edf.1"%(scanname,counter,idet,scannumber,f),\
                    "%s_%s%02d_%04d_%s.edf"%(scanname,counter,idet,scannumber,f),\
                    "%s_%s%02d_%04d_%s.edf.1"%(scanname,counter,idet,scannumber,f)]
    else:
        filename = ["%s_%s_%04d_%s.edf"%(scanname,counter,scannumber,f),\
                    "%s_PUZ_%s_%04d_%s.edf"%(scanname,counter,scannumber,f),\
                    "%s_%s_%04d_%s.edf.1"%(scanname,counter,scannumber,f)]

    n = 0
    for name in filename: # loop over possibilities
        counterfile = os.path.join(sourcepath,name)
        files = glob(counterfile)
        files = [f for f in files if os.stat(f).st_size > 0]
        n = len(files)
        if n!=0:
            break
        
    if getcount:
        return n
    else:
        if n==0:
            return None
        else:
            return counterfile

def detectorname(config,ndet,idet,counter=False):
    if ndet == 0:
        return ""
    elif config["addbeforefitting"] and ndet > 1 and not counter:
        return "detectorsum"
    else:
        return "detector%d"%idet

def dimensions(config):
    stackdim = config["stackdim"]
    if stackdim == 0:
        imgdim = [1,2]
    elif stackdim == 1:
        imgdim = [0,2]
    else:
        imgdim = [0,1]
    return stackdim,imgdim

def getscanpositions(config,header):
    """Get scan dimensions from header
    """
    result = {"name":"unknown"}
    
    if "scanlabel" in config:
        if config["scanlabel"] in header:
            o = spec.cmd_parser()
            result = o.parsezapimage(header[config["scanlabel"]])
    elif "fastlabel" in config and "slowlabel" in config:
        o = spec.edfheader_parser(fastlabel=config["fastlabel"],slowlabel=config["slowlabel"])
        result = o.parsezapimage(header)

    if result["name"]=="zapimage":
        sfast = {"name":result["motfast"],"data":spec.zapline_values(result["startfast"],result["endfast"],result["npixelsfast"])} 
        sslow = {"name":result["motslow"],"data":spec.ascan_values(result["startslow"],result["endslow"],result["nstepsslow"])} 
        return (sfast["name"],sslow["name"],sfast,sslow)
    else:
        raise RuntimeError("Scan command cannot be parsed for motor positions.")

def getimagestacks(config):
    """Get image stacks (counters, ROI's, fitted maps)

    Args:

    Returns:
        tuple

        The first element contains the image stack:
            stacks = {"counters":{"name1":lstack1,"name2":lstack2,...},
                      "det0":{"name3":lstack3,"name4":lstack4,...},
                      "det1":{"name3":lstack3,"name4":lstack4,...},...}
            lstack: an image stack given as a list of strings (filenames)

        The second element is a list with three elements which contains
        the axis values of the stack:
            stackaxes = [{"name":"name1","data":np.array},
                         {"name":"name2","data":np.array},
                         {"name":"name3","data":np.array}]

        The third element is a dictionary of static coordinates:
            coordinates = {"varname1":np.array, "varname2":np.array, ...}
    """

    

    # Prepare data
    npaths = len(config["sourcepath"])
    if npaths != len(config["scanname"]):
        raise ValueError("Number of scan names must be the same as number of source paths.")
    if npaths != len(config["scannumbers"]):
        raise ValueError("Number of scan numbers must be the same as number of source paths.")

    stacks = {}
    stackaxes = [None]*3
    stackdim,imgdim = dimensions(config)
    coordinates = {}
    nscanstot = 0
    iscanoffset = 0
    for ipath in range(npaths):
        nscanstot += len(config["scannumbers"][ipath])

    # DT correction, sum detector and fit
    nfilestofit = None
    for ipath in range(npaths):
        counterpath = os.path.abspath(os.path.join(config["sourcepath"][ipath],config["counter_reldir"]))
        xrfpath = config["sourcepath"][ipath]
        
        scanname = config["scanname"][ipath]
        nscans = len(config["scannumbers"][ipath])

        for iscan in range(nscans):
            scannumber = config["scannumbers"][ipath][iscan]
            
            logger.info("{} ({})...".format(scanname,scannumber))

            # Extract stack axes values
            stackvalue = np.nan
            for metacounter in config["metacounters"]:
                try:
                    if "xmap" in metacounter:
                        idet = 0
                    else:
                        idet = None
                    metafilename = filecounter(counterpath,scanname,metacounter,scannumber,idet=idet)
                    metafile = fabio.open(metafilename)
                    header = metafile.header

                    if iscan == 0 and ipath == 0:
                        try:
                            motfast,motslow,sfast,sslow = getscanpositions(config,header)
                        except:
                            motfast = "fast"
                            motslow = "slow"

                            tmp = int(header["Dim_2"])
                            sfast = {"name":motfast,"data":np.arange(tmp)}

                            tmp = filecounter(counterpath,scanname,metacounter,scannumber,idet=idet,getcount=True)
                            sslow = {"name":motslow,"data":np.arange(tmp)}
                            
                        stackaxes[imgdim[1]] = sfast
                        stackaxes[imgdim[0]] = sslow
                        stackaxes[stackdim] = {"name":str(config["stacklabel"]),"data":np.full(nscanstot,np.nan,dtype=np.float32)}

                        coordinates = {}
                        for mot in config["coordinates"]:
                            if mot != motfast and mot != motslow and mot in header:
                                coordinates[mot] = np.full(nscanstot,np.nan)
                    
                    # Get coordinates
                    for mot in coordinates:
                        if mot in header:
                            coordinates[mot][iscanoffset+iscan] = np.float(header[mot])
                            
                    # Get stack value
                    if config["stacklabel"] in header:
                        stackvalue = np.float(header[config["stacklabel"]])
                    break
                except:
                    if metafilename is None:
                        metafilename = filecounter_nocheck(counterpath,scanname,metacounter,scannumber,idet=idet)
                    logger.exception("Something wrong with extracting info from meta file {}.".format(metafilename))

            # Stack axis
            if stackaxes[stackdim] is None:
                raise IOError("Metacounter files are not present, corrupted or not the right format.")
            stackaxes[stackdim]["data"][iscanoffset+iscan] = stackvalue

            # DT correction and sum detector (if more than 1)
            if config["dtcor"]:
                parsename = "dtcor"
            else:
                parsename = "raw"
            parsename = "%%0%dd_%s"%(np.int(np.floor(np.log10(npaths)))+1,parsename)%(ipath)

            filestofit,detnums = parse_xia_esrf(xrfpath,scanname,
                    scannumber,config["outdatapath"],parsename,add=config["addbeforefitting"],
                    exclude_detectors=config["exclude_detectors"],deadtime=config["dtcor"],
                    willbefitted=config["fit"],deadtimeifsingle=config["dtcorifsingle"])
            ndet = len(detnums)
            noxia = ndet==0

            # XIA files may be missing but maybe the xia counters are there
            if noxia:
                detcounters = [ctr for ctr in config["counters"] if "xmap" in ctr]
                if len(detcounters)!=0:
                    idet = 0
                    while True:                      
                        tmp = filecounter(counterpath,scanname,detcounters[0],scannumber,idet=idet,getcount=True)
                        if tmp==0:
                            break
                        idet += 1
                    ndet = idet
                    detnums = range(ndet)

            # Allocate first level of the stacks
            if iscan == 0 and ipath == 0:
                counters = [ctr for ctr in config["counters"] if "xmap" not in ctr]
                detcounters = [ctr for ctr in config["counters"] if "xmap" in ctr]

                for idet in detnums:
                    stacks[detectorname(config,ndet,idet)] = {}
                    if len(detcounters)!=0:
                        stacks[detectorname(config,ndet,idet,counter=True)] = {}

                stacks["counters"] = {}
                for counter in counters:
                    stacks["counters"][counter] = [""]*nscanstot

            # Loop over the detectors
            for i in range(ndet):
                idet = detnums[i]

                # Fit data
                if not config["addbeforefitting"] or i==0:
                    if config["fit"] and not noxia:
                        if len(filestofit[i])!= 0:
                            if len(config["detectorcfg"])==1:
                                cfg = config["detectorcfg"][0]
                            else:
                                cfg = config["detectorcfg"][i]

                            # Fit spectra resulting in images
                            if config["addbeforefitting"]:
                                outname = "%s_%s_xiaS1_%04d_0000"%(scanname,parsename,scannumber)
                            else:
                                outname = "%s_%s_xia%02d_%04d_0000"%(scanname,parsename,idet,scannumber)
                            files, labels = fitter(filestofit[i],
                                                   config["outfitpath"],outname,cfg,stackvalue,
                                                   fast=config["fastfitting"],mlines=config["mlines"])

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
                    stacks[detname][counter][iscanoffset+iscan] = filecounter(counterpath,scanname,counter,scannumber,idet=idet)

            # Append counters
            for counter in counters:
                stacks["counters"][counter][iscanoffset+iscan] = filecounter(counterpath,scanname,counter,scannumber)
        
        iscanoffset += nscans

    # Sort stack on stack axis value
    ind = np.argsort(stackaxes[stackdim]["data"],kind='mergesort')
    stackaxes[stackdim]["data"] = stackaxes[stackdim]["data"][ind]
    for mot in coordinates:
        coordinates[mot] = coordinates[mot][ind]
    for s in stacks:
        group = stacks[s]
        for lstack in group:
            group[lstack] = [group[lstack][i] for i in ind]

    return stacks,stackaxes,coordinates

def exportgroups(f,stacks,keys,axes,stackdim,imgdim,sumgroups=False):
    """Export groups of EDF stacks, summated or not
    """

    logger = logging.getLogger(__name__)

    if sumgroups:
        sumname = "detectorsum"
        if sumname in f:
            grpsum = f[sumname]
        else:
            grpsum = nexus.newNXentry(f,sumname)
        if sumname not in stacks:
            stacks[sumname] = {}

    dim = [0,0,0]

    for k1 in keys: # detector or counter group
        bgrpissum = k1 is not "counters" and sumgroups
        if bgrpissum:
            grp = grpsum
        else:
            if k1 in f:
                grp = f[k1]
            else:
                grp = nexus.newNXentry(f,k1)

        for k2 in stacks[k1].keys(): # stack subgroup (Al-K, S-K, xmap_x1c, ...)
            nxdatagrp = None

            # skip when already an Nxdata group (previous call to exportgroups)
            if not hasattr(stacks[k1][k2],"__iter__"):
                continue

            # loop over the stack images
            # stacks[k1][k2]: list of files
            nscans = len(stacks[k1][k2])
            for iscan in range(nscans):
                filename = stacks[k1][k2][iscan]
                if filename is None:
                    continue
                data = fabio.open(filename).data

                if k2 in grp:
                    dset = grp[k2][grp[k2].attrs["signal"]]
                else:
                    nxdatagrp = nexus.newNXdata(grp,k2,"")
                    
                    # Allocate dataset
                    if dim==[0,0,0]:
                        dim[imgdim[0]] = data.shape[0]
                        dim[imgdim[1]] = data.shape[1]
                        dim[stackdim] = nscans
                    dset = nexus.createNXdataSignal(nxdatagrp,shape=dim,chunks = True,dtype = np.float32)
                    if bgrpissum:
                        dset[:] = 0
                    else:
                        dset[:] = np.nan

                    # Link axes to the new NXdata group
                    nexus.linkaxes(f,axes,[nxdatagrp])

                # Some rows too much or rows missing:
                if dim[imgdim[0]] > data.shape[0]:
                    data = np.pad(data,((0,dim[imgdim[0]]-data.shape[0]),(0,0)),'constant',constant_values=0)
                elif dim[imgdim[0]] < data.shape[0]:
                    data = data[0:dim[imgdim[0]],:]

                if stackdim == 0:
                    if bgrpissum:
                        dset[iscan,...] += data
                    else:
                        dset[iscan,...] = data
                elif stackdim == 1:
                    if bgrpissum:
                        dset[:,iscan,:] += data
                    else:
                        dset[:,iscan,:] = data
                else:
                    if bgrpissum:
                        dset[...,iscan] += data
                    else:
                        dset[...,iscan] = data

            # Add NX entry to stacks
            if nxdatagrp is None:
                del stacks[k1][k2]
            else:
                if bgrpissum:
                    stacks[sumname][k2] = nxdatagrp.name
                else:
                    stacks[k1][k2] = nxdatagrp.name
            
        if bgrpissum:
            del stacks[k1]

    return stacks

def exportimagestacks(config,stacks,stackaxes,coordinates,jsonfile):
    """Export EDF stacks to HDF5
    """
    f = nexus.File(config["hdf5output"],mode='w')

    # Save stack axes values
    axes = nexus.createaxes(f,stackaxes)

    stackdim,imgdim = dimensions(config)

    # Save sum of detectors if requested
    keys = [k for k in stacks.keys() if "detector" in k and "sum" not in k]
    if config["addafterfitting"] and len(keys)>1:
        stacks = exportgroups(f,stacks,keys,axes,stackdim,imgdim,sumgroups=True)

    # Save other groups (not included in the sum)
    if len(stacks) > 0:
        stacks = exportgroups(f,stacks,stacks.keys(),axes,stackdim,imgdim)

    # Save coordinates
    coordgrp = nexus.newNXentry(f,"coordinates")
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

    Returns:
        tuple

        The first element contains the image stack:
            stacks = {"counters":{"name1":lstack1,"name2":lstack2,...},
                      "det0":{"name3":lstack3,"name4":lstack4,...},
                      "det1":{"name3":lstack5,"name4":lstack6,...},...}
            lstack: an image stack given as an NXdata path

        The second element is a list with three elements which contains
        the axis of the stack:
            axes = [{"name":"name1","fullname":"/axes/name1/data"},
                    {"name":"name2","fullname":"/axes/name2/data"},
                    {"name":"name3","fullname":"/axes/name3/data"}]
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
        create_hdf5_imagestacks(sys.argv[1])



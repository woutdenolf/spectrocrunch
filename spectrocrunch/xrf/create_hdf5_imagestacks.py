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

import json
import numpy as np
import logging

from ..io import xiaedf
from ..io import edf
from ..io import spec
from ..io import nexus
from ..xrf.fit import PerformBatchFit as fitter

logger = logging.getLogger(__name__)

def axesindices(config):
    """Get stack and image axes indices
    
    Args:
        config(dict):
        
    Returns:
        stackdim(num)
        imgdim(list)
    """
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
            result = o.parse(header[config["scanlabel"]])
    elif "fastlabel" in config and "slowlabel" in config:
        o = spec.edfheader_parser(fastlabel=config["fastlabel"],slowlabel=config["slowlabel"])
        result = o.parse(header)

    if result["name"]=="zapimage":
        sfast = {"name":result["motfast"],"data":spec.zapline_values(result["startfast"],result["endfast"],result["npixelsfast"])} 
        sslow = {"name":result["motslow"],"data":spec.ascan_values(result["startslow"],result["endslow"],result["nstepsslow"])} 
    else:
        sfast = {"name":"fast","data":np.arange(header["Dim_2"])}
        sslow = {"name":"slow","data":np.arange(header["Dim_1"])}
 
    return (sfast["name"],sslow["name"],sfast,sslow)

def detectorname(detector):
    if detector:
        if detector.isdigit():
            name = "detector{:d}".format(int(detector))
        elif detector.startswith("S"):
            name = "detectorsum"
        else:
            raise "Unexpected detector name {}".format(detector)
    else:
        name = "counter"
    return name

def getimagestacks(config):
    """Get image stacks (counters, ROI's, fitted maps)

    Args:

    Returns:
        stacks(dict): {"counters":{"name1":filenames1,"name2":filenames2,...},
                       "det0":{"name3":filenames3,"name4":filenames4,...},
                       "det1":{"name3":filenames3,"name4":filenames4,...},
                       ...}
        stackaxes(dict): [{"name":"name1","data":np.array},
                          {"name":"name2","data":np.array},
                          {"name":"name3","data":np.array}]

        coordinates(dict): {"motorname1":np.array, "motorname2":np.array, ...}
    """

    # Check data
    npaths = len(config["sourcepath"])
    if npaths != len(config["scanname"]):
        raise ValueError("Number of scan names must be the same as number of source paths.")
    if npaths != len(config["scannumbers"]):
        raise ValueError("Number of scan numbers must be the same as number of source paths.")

    # Initialize result
    stacks = {}
    stackaxes = [None]*3
    coordinates = {}
    
    stackdim,imgdim = axesindices(config)
    
    # Get image stack
    xiastackraw = xiaedf.xiastack_mapnumbers(config["sourcepath"],config["scanname"],config["scannumbers"])
    
    # Exclude detectors
    xiastackraw.skipdetectors(config["exclude_detectors"])
    
    # Counter directory relative to the XIA files
    xiastackraw.counter_reldir(config["counter_reldir"])
    
    # Stack shape
    nstack, nrow, ncol, nchan, ndet = xiastackraw.dshape

    # Deadtime correction
    if config["dtcor"]:
        if ndet==1:
            dtcor = config["dtcorifsingle"]
        else:
            dtcor = config["fit"]
    else:
        dtcor = False
    xiastackraw.dtcor(dtcor)
    
    # Add detectors
    adddet = config["fit"] and config["addbeforefitting"] and ndet>1
    xiastackraw.detectorsum(adddet)
    
    # Create new spectra when needed
    if dtcor or adddet or True: # TODO: remove True
        if dtcor:
            label = "dtcor"
        else:
            label = "raw"
        radix = ["{}_{}".format(radix,label) for radix in config["scanname"]]
        xiastackproc = xiaedf.xiastack_mapnumbers(config["outdatapath"],radix,config["scannumbers"])
        xiastackproc.overwrite(True)
        if adddet:
            xialabels = ["xiaS1"]
        else:
            xialabels = ["xia{:02d}".format(det) for det in xiastackraw.xiadetectorselect_numbers(range(ndet))]
        
        xiastackproc.save(xiastackraw.data,xialabels)
        nstack, nrow, ncol, nchan, ndet = xiastackproc.dshape
    else:
        xiastackproc = xiastackraw
    
    # Check counters
    counters = set(xiastackraw.counterbasenames())
    metacounters = counters.intersection(config["metacounters"])
    counters = counters.intersection(config["counters"])
    
    if metacounters:
        metacounters = next(iter(metacounters))
    elif "xia" in config["metacounters"]:
        metacounters = "xia"
    else:
        logger.exception("Metacounters for {} are not found.".format(xiastackraw)) 

    # Extract metadata and counters from raw stack
    for imageindex,xiaimage in enumerate(xiastackraw):
        binit = imageindex==0

        if metacounters=="xia":
            files = xiaimage.statfilenames()
        else:
            files = xiaimage.ctrfilenames(metacounters)

        if files:
            header = edf.edfimage(files[0]).header
            motfast,motslow,sfast,sslow = getscanpositions(config,header)
            
            # Prepare axes and coordinates
            if binit:
                for mot in config["coordinates"]:
                    if mot != motfast and mot != motslow and mot in header:
                        coordinates[mot] = np.full(nstack,np.nan)
                        
                stackaxes[imgdim[1]] = sfast
                stackaxes[imgdim[0]] = sslow
                stackaxes[stackdim] = {"name":str(config["stacklabel"]),"data":np.full(nstack,np.nan,dtype=np.float32)}

            # Add coordinates
            for mot in coordinates:
                if mot in header:
                    coordinates[mot][imageindex] = np.float(header[mot])
                              
            # Add stack value
            if config["stacklabel"] in header:
                stackaxes[stackdim]["data"][imageindex] = np.float(header[config["stacklabel"]])
        else:
            logger.exception("Metacounters for {} are not found.".format(xiaimage)) 
    
        # Counters
        files = xiaimage.ctrfilenames(counters)
        files = xiaedf.xiagroupdetectors(files)
        for detector,v1 in files.items():
            name = detectorname(detector)
            
            # Prepare list of files
            if binit:
                stacks[name] = {}
                for ctr in v1:
                    stacks[name][ctr] = [""]*nstack

            for ctr,f in v1.items():
                # Add counter file
                stacks[name][ctr][imageindex] = f[0]

    # Fit data and add elemental maps
    if config["fit"]:
        if len(config["detectorcfg"])==1:
            fitcfg = config["detectorcfg"]*ndet
        else:
            fitcfg = config["detectorcfg"]

        for imageindex,xiaimage in enumerate(xiastackproc):
            binit = imageindex==0
                
            filestofit = xiaimage.datafilenames()
            filestofit = xiaedf.xiagroupdetectors(filestofit)
            for detector,cfg in zip(filestofit,fitcfg):
                # Fit
                outname = "{}_xia{}_{:04d}_0000".format(xiaimage.radix,detector,xiaimage.mapnum)
                energy = stackaxes[stackdim]["data"][imageindex]
                files, labels = fitter(filestofit[detector]["xia"],
                                       config["outfitpath"],outname,cfg,energy,
                                       fast=config["fastfitting"],mlines=config["mlines"])
                
                # Prepare list of files    
                name = detectorname(detector)
                if binit:
                    if name not in stacks:
                        stacks[name] = {}
                    for label in labels:
                        stacks[name][label] = [""]*nstack
                 
                # Add file name           
                for f,label in zip(files,labels):
                    stacks[name][label][imageindex] = f
                    
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
                data = edf.edfimage(filename).data

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

    stackdim,imgdim = axesindices(config)

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
        

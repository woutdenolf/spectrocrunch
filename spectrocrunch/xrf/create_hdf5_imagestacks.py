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
import re
import shutil
import os
import collections

from ..io import xiaedf
from ..io import edf
from ..io import spec
from ..io import nexus
from ..xrf.fit import PerformBatchFit
from ..common import units

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

def getscanparameters(config,header):
    """Get scan dimensions from header
    """
    result = {"name":"unknown"}
    
    if "speccmdlabel" in config:
        if config["speccmdlabel"] in header:
            o = spec.cmd_parser()
            result = o.parse(header[config["speccmdlabel"]])
    elif "fastlabel" in config and "slowlabel" in config:
        o = spec.edfheader_parser(fastlabel=config["fastlabel"],slowlabel=config["slowlabel"])
        result = o.parse(header)

    if result["name"]=="zapimage":
        sfast = {"name":result["motfast"],"data":spec.zapline_values(result["startfast"],result["endfast"],result["npixelsfast"])} 
        sslow = {"name":result["motslow"],"data":spec.ascan_values(result["startslow"],result["endslow"],result["nstepsslow"])} 
        if "time" in result:
            time = units.umagnitude(result["time"],units="s")
        else:
            time = np.nan
    else:
        logger.warning("No motor positions in header (using pixels).")
        sfast = {"name":"fast","data":np.arange(int(header["Dim_2"]))}
        sslow = {"name":"slow","data":np.arange(int(header["Dim_1"]))}
        time = np.nan
        
    return sfast["name"],sslow["name"],sfast,sslow,time

def detectorname(detector):
    # detector: "00", "S0"
    if detector:
        if detector.isdigit():
            name = "detector{:d}".format(int(detector))
        elif detector.startswith("S"):
            name = "detectorsum"
        else:
            raise "Unexpected detector name {}".format(detector)
    else:
        name = "counters"
    return name

def createimagestacks(config,qxrfgeometry=None):
    """Get image stacks (counters, ROI's, fitted maps)

    Args:
        config(dict):
        qxrfgeometry(Optional(object)):
        
    Returns:
        stacks(dict): {"counters":{"name1":filenames1,"name2":filenames2,...},
        "det0":{"name3":filenames3,"name4":filenames4,...},
        "det1":{"name3":filenames3,"name4":filenames4,...},
        ...}
        
        stackaxes(dict): [{"name":"name1","data":np.array},
        {"name":"name2","data":np.array},
        {"name":"name3","data":np.array}]

        stackinfo(dict): {"motorname1":np.array, "motorname2":np.array, ...}
        
        procinfo(dict): processing info
    """
    
    procinfo = {}
    
    # Check data
    npaths = len(config["sourcepath"])
    if npaths != len(config["scanname"]):
        raise ValueError("Number of scan names must be the same as number of source paths.")
    if npaths != len(config["scannumbers"]):
        raise ValueError("Number of scan numbers must be the same as number of source paths.")

    # Initialize result
    stacks = {}
    stackaxes = [None]*3
    stackinfo = {}
    
    stackdim,imgdim = axesindices(config)
    
    # Get image stack
    xiastackraw = xiaedf.xiastack_mapnumbers(config["sourcepath"],config["scanname"],config["scannumbers"])
    if xiastackraw.isempty:
        raise IOError("Cannot find data: {}".format(xiastackraw.filedescription))
    nstack, nrow, ncol, nchan, ndetorg = xiastackraw.dshape

    # Exclude detectors
    xiastackraw.skipdetectors(config["exclude_detectors"])
    xiastackraw.keepdetectors(config["include_detectors"])
    nstack, nrow, ncol, nchan, ndet = xiastackraw.dshape
    
    # Counter directory relative to the XIA files
    xiastackraw.counter_reldir(config["counter_reldir"])
    
    # Processing
    addbefore = config["addbeforefitting"] and ndet>1
    addafter = config["addafterfitting"] and ndet>1 and not addbefore
    fluxnorm = qxrfgeometry is not None
    dtcor = config["dtcor"] and (config["dtcorifsingle"] or ndet>1 or fluxnorm)
    fit = config["fit"]
    
    procinfo["dtneeded"] = config["dtcor"] and not dtcor
    procinfo["fluxnorm"] = fluxnorm
    
    xiastackraw.detectorsum(addbefore)
    xiastackraw.dtcor(dtcor)

    # Check counters
    countersfound = set(xiastackraw.counterbasenames())
    counters = countersfound.intersection(config["counters"])
    if config["metadata"]=="xia":
        metacounters = "xia"
    else:
        if countersfound:
            metacounters = next(iter(countersfound))
        else:
            logger.warning("Metacounters for {} are not found".format(xiastackraw)) 
            metacounters = []
            
    # Extract metadata and counters from raw stack
    for imageindex,xiaimage in enumerate(xiastackraw):
        binit = imageindex==0

        if metacounters=="xia":
            files = xiaimage.statfilenames()
        else:
            files = xiaimage.ctrfilenames(metacounters)
        if files:
            header = edf.edfimage(files[0]).header
        else:
            logger.warning("Metacounters for {} are not found".format(xiaimage)) 
            header = {"Dim_1":ncol, "Dim_2":nrow}
        motfast,motslow,sfast,sslow,time = getscanparameters(config,header)

        # Prepare axes and stackinfo
        if binit:
            for mot in config["stackinfo"]:
                if mot != motfast and mot != motslow and mot in header:
                    stackinfo[mot] = np.full(nstack,np.nan)
                    
            stackaxes[imgdim[1]] = sfast
            stackaxes[imgdim[0]] = sslow
            stackaxes[stackdim] = {"name":str(config["stacklabel"]),"data":np.full(nstack,np.nan,dtype=np.float32)}
            stackinfo["expotime"] = np.full(nstack,np.nan,dtype=np.float32)
            stackinfo["sampledetdistance"] = np.full(nstack,np.nan,dtype=np.float32)

        # Add stackinfo
        for mot in stackinfo:
            if mot in header:
                stackinfo[mot][imageindex] = np.float(header[mot])
                          
        # Add stack value
        if config["stacklabel"] in header:
            stackaxes[stackdim]["data"][imageindex] = np.float(header[config["stacklabel"]])
        else:
            logger.warning("No stack counter in header (set to NaN)")
        
        # Add time
        stackinfo["expotime"][imageindex] = time
        
        # Counters
        files = xiaimage.ctrfilenames(counters)
        files = xiaedf.xiagroupdetectors(files)
        for detector,v1 in files.items():
            name = detectorname(detector)
            
            # Prepare list of files
            if binit:
                stacks[name] = collections.OrderedDict()
                for ctr in v1:
                    stacks[name][ctr] = [""]*nstack
            for ctr,f in v1.items():
                # Add counter file
                stacks[name][ctr][imageindex] = f[0]
                
    # Create new spectra when needed
    if dtcor or addbefore or fluxnorm:
    
        # Set normalizer for each image separately
        if fluxnorm:
            stackinfo["refflux"] = np.full(nstack,np.nan,dtype=np.float32)
            stackinfo["refexpotime"] = np.full(nstack,np.nan,dtype=np.float32)
            stackinfo["activearea"] = np.full(nstack,qxrfgeometry.xrfgeometry.detector.activearea)
            stackinfo["anglein"] = np.full(nstack,qxrfgeometry.xrfgeometry.anglein)
            stackinfo["angleout"] = np.full(nstack,qxrfgeometry.xrfgeometry.angleout)
            
            for imageindex,xiaimage in enumerate(xiastackraw):
                energy = stackaxes[stackdim]["data"][imageindex]
                if not np.isnan(energy):
                    time = stackinfo["expotime"][imageindex]
                    if np.isnan(time):
                        time = None
                        
                    xrfnormop,\
                    stackinfo["refflux"][imageindex],\
                    stackinfo["refexpotime"][imageindex],\
                    stackinfo["expotime"][imageindex]\
                     =qxrfgeometry.xrfnormop(energy,expotime=time)
                    
                    xiaimage.localnorm(config["fluxcounter"],func=xrfnormop)

                    pos = stackinfo["sampledetdistance"][imageindex]
                    if not np.isnan(pos):
                        qxrfgeometry.setxrfposition(pos)
                    stackinfo["sampledetdistance"][imageindex] = qxrfgeometry.getxrfdistance()
        
        label = ""
        if dtcor:
            label = label+"dt"
        if fluxnorm:
            label = label+"fl"
        if label:
            label = label+"cor"
            radix = ["{}_{}".format(radix,label) for radix in config["scanname"]]
        else:
            radix = config["scanname"]
        
        # not necesarry but clean in case of re-runs
        if os.path.isdir(config["outdatapath"]):
            shutil.rmtree(config["outdatapath"])
        xiastackproc = xiaedf.xiastack_mapnumbers(config["outdatapath"],radix,config["scannumbers"])
        xiastackproc.overwrite(True)
        
        if addbefore:
            xialabels = ["xiaS1"]
        else:
            xialabels = xiastackraw.xialabels_used

        logger.info("Creating corrected XRF spectra ...")
        #xiastackproc.save(xiastackraw.data,xialabels=xialabels) # all in memory
        xiastackproc.save(xiastackraw,xialabels=xialabels) # least memory usage possible
        nstack, nrow, ncol, nchan, ndet = xiastackproc.dshape
    else:
        xiastackproc = xiastackraw
    
    # I0/It stacks
    if fluxnorm:
        energy = stackaxes[stackdim]["data"][imageindex]
        if not np.isnan(energy) and ("fluxcounter" in config or "transmissioncounter" in config):
            for imageindex in range(nstack):
                name = detectorname(None)
                time = stackinfo["refexpotime"][imageindex]
                if "fluxcounter" in config:
                    op,_ = qxrfgeometry.I0op(energy,expotime=time)
                    if "calc_flux0" not in stacks[name]:
                        stacks[name]["calc_flux0"] = [""]*nstack
                    stacks[name]["calc_flux0"][imageindex] = {"args":[(name,config["fluxcounter"])],"func":op}
                if "transmissioncounter" in config:
                    op,_ = qxrfgeometry.Itop(energy,expotime=time)
                    if "calc_fluxt" not in stacks[name]:
                        stacks[name]["calc_fluxt"] = [""]*nstack
                        stacks[name]["calc_transmission"] = [""]*nstack
                    stacks[name]["calc_fluxt"][imageindex] = {"args":[(name,config["transmissioncounter"])],"func":op}
                    stacks[name]["calc_transmission"][imageindex] = {"args":[(name,"calc_fluxt"),(name,"calc_flux0")],"func":lambda a,b: -np.log(a/b.astype(np.float32))}
                    
    # Fit data and add elemental maps
    if fit:
        logger.info("Fit XRF spectra ...")
        
        # not necesarry but clean in case of re-runs
        if os.path.isdir(config["outfitpath"]):
            shutil.rmtree(config["outfitpath"])
                    
        if len(config["detectorcfg"])==1:
            fitcfg = config["detectorcfg"]*ndet
        else:
            fitcfg = config["detectorcfg"]
            if len(fitcfg)!=ndet:
                raise RuntimeError("You need {} configuration files, {} provides.".format(ndet,len(fitcfg)))
                
        for imageindex,xiaimage in enumerate(xiastackproc):
            binit = imageindex==0
            
            if fluxnorm:
                quant = {"time":stackinfo["refexpotime"][imageindex],\
                        "flux":stackinfo["refflux"][imageindex],\
                        "area":stackinfo["activearea"][imageindex],\
                        "anglein":stackinfo["anglein"][imageindex],\
                        "angleout":stackinfo["angleout"][imageindex],\
                        "distance":stackinfo["sampledetdistance"][imageindex]}
            else:
                quant = {}
            
            filestofit = xiaimage.datafilenames_used()
            filestofit = xiaedf.xiagroupdetectors(filestofit)
            for detector,cfg in zip(filestofit,fitcfg):
                # Fit
                outname = "{}_xia{}_{:04d}_0000".format(xiaimage.radix,detector,xiaimage.mapnum)
                energy = stackaxes[stackdim]["data"][imageindex]

                files, labels = PerformBatchFit(filestofit[detector]["xia"],
                                       config["outfitpath"],outname,cfg,energy,
                                       fast=config["fastfitting"],mlines=config["mlines"],quant=quant)
                
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
    
    # Detector sum
    if addafter:
        detectors = [k for k in stacks.keys() if re.match("^detector([0-9]+)$",k)]
        if "detectorsum" not in stacks:
            stacks["detectorsum"] = {}
        
        for k1 in detectors:
            for k2 in stacks[k1]:
                if k2 not in stacks["detectorsum"]:
                    stacks["detectorsum"][k2] = [[]]*nstack
                    
                for imageindex in range(nstack):
                    stacks["detectorsum"][k2][imageindex].append(stacks[k1][k2][imageindex])
            stacks.pop(k1)
        
    # Sort stack on stack axis value
    ind = np.argsort(stackaxes[stackdim]["data"],kind='mergesort')
    stackaxes[stackdim]["data"] = stackaxes[stackdim]["data"][ind]
    for mot in stackinfo:
        stackinfo[mot] = stackinfo[mot][ind]
    for k1 in stacks:
        group = stacks[k1]
        for k2 in group:
            group[k2] = [group[k2][i] for i in ind]

    return stacks,stackaxes,stackinfo,procinfo

def exportgroups(f,stacks,axes,stackdim,imgdim,stackshape,proc):
    """Export groups of EDF stacks, summated or not
    """

    for k1 in stacks: # detector or counter group
        if k1 in f:
            grp = f[k1]
        else:
            grp = nexus.newNXentry(f,k1)

        for k2 in stacks[k1]: # stack subgroup (Al-K, S-K, xmap_x1c, ...)
            # skip when already an Nxdata group (previous call to exportgroups)
            if not hasattr(stacks[k1][k2],"__iter__"):
                continue

            # remove empty stack
            nscans = len(stacks[k1][k2])
            if all(stacks[k1][k2][iscan] is None for iscan in range(nscans)):
                stacks[k1].pop(k2)
                continue

            # loop over the stack images
            # stacks[k1][k2]: list of files or process dictionaries
            nxdatagrp = None
            for iscan in range(nscans):
                datainfo = stacks[k1][k2][iscan]
                
                # Get data (1 image from the stack)
                if proc=="calc":
                    if not isinstance(datainfo,dict):
                        break
                    datainfo = stacks[k1][k2][iscan]

                    args = []
                    for k1b,k2b in datainfo["args"]:
                        grp2 = f[k1b]
                        dset = grp2[k2b][grp2[k2b].attrs["signal"]]
                        if stackdim == 0:
                            data = dset[iscan,...]
                        elif stackdim == 1:
                            data = dset[:,iscan,:]
                        else:
                            data = dset[...,iscan]
                        args.append(data)

                    data = datainfo["func"](*args)
                else:
                    if isinstance(datainfo,dict):
                        break
                    if not isinstance(datainfo,list):
                        datainfo = [datainfo]
                    data = sum(edf.edfimage(filename).data for filename in datainfo)

                # Get destination for data
                if k2 in grp:
                    dset = grp[k2][grp[k2].attrs["signal"]]
                else:
                    logger.info("saving {}/{}".format(k1,k2))
                    nxdatagrp = nexus.newNXdata(grp,k2,"")
                    
                    # Allocate dataset
                    if stackshape==[0,0,0]:
                        stackshape[imgdim[0]] = data.shape[0]
                        stackshape[imgdim[1]] = data.shape[1]
                        stackshape[stackdim] = nscans
                    dset = nexus.createNXdataSignal(nxdatagrp,shape=stackshape,chunks = True,dtype = np.float32)
                    dset[:] = np.nan

                    # Link axes to the new NXdata group
                    nexus.linkaxes(f,axes,[nxdatagrp])

                # Some rows too much or rows missing:
                if stackshape[imgdim[0]] > data.shape[0]:
                    data = np.pad(data,((0,stackshape[imgdim[0]]-data.shape[0]),(0,0)),'constant',constant_values=0)
                elif stackshape[imgdim[0]] < data.shape[0]:
                    data = data[0:stackshape[imgdim[0]],:]

                # Save data
                if stackdim == 0:
                    dset[iscan,...] = data
                elif stackdim == 1:
                    dset[:,iscan,:] = data
                else:
                    dset[...,iscan] = data

            else:
                if nxdatagrp is not None:
                    # Replace subgroup k2 filename or calcinfo with NXentry name in stack
                    stacks[k1][k2] = nxdatagrp.name

    return stacks,stackshape

def exportimagestacks(config,stacks,stackaxes,stackinfo,procinfo):
    """Export EDF stacks to HDF5
    """
    
    with nexus.File(config["hdf5output"],mode='w') as f:
        # Save stack axes values
        axes = nexus.createaxes(f,stackaxes)

        # Save groups
        stackdim,imgdim = axesindices(config)
        stackshape = [0,0,0]
        stacks,stackshape = exportgroups(f,stacks,axes,stackdim,imgdim,stackshape,"raw")
        stacks,stackshape = exportgroups(f,stacks,axes,stackdim,imgdim,stackshape,"calc")

        # Save stackinfo
        coordgrp = nexus.newNXentry(f,"stackinfo")
        for k in stackinfo:
            coordgrp[k] = stackinfo[k]

        # Add processing info
        nexus.addinfogroup(f,"fromraw",procinfo)

    logger.info("Saved {}".format(config["hdf5output"]))

    return axes
    
def create_hdf5_imagestacks(jsonfile,qxrfgeometry=None):
    """Convert scanning data (XIA spectra + counters) to an HDF5 file:
        groups which contain NXdata classes
        3 axes datasets on the main level

    Returns:
        stacks(dict): {"counters": {"name1":lstack1,"name2":lstack2,...},
        "detector0":{"name3":lstack3,"name4":lstack4,...},
        "detector1":{"name3":lstack5,"name4":lstack6,...},...}
        lstack: an image stack given as an NXdata path

        axes(list): [{"name":"name1","fullname":"/axes/name1/data"},
        {"name":"name2","fullname":"/axes/name2/data"},
        {"name":"name3","fullname":"/axes/name3/data"}]
    """
    # Processing configuration
    with open(jsonfile,'r') as f:
        config = json.load(f)

    # Raw or pre-processed data (e.g. DT correction, fitting)
    stacks,stackaxes,stackinfo,procinfo = createimagestacks(config,qxrfgeometry=qxrfgeometry)

    # Export EDF stacks to HDF5 stacks
    procinfo["config"] = jsonfile
    axes = exportimagestacks(config,stacks,stackaxes,stackinfo,procinfo)

    return stacks,axes,[procinfo]

if __name__ == '__main__':
    import sys
    if len(sys.argv)>=2:
        create_hdf5_imagestacks(sys.argv[1])
        

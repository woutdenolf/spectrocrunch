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
import h5py

from ..xrf.create_hdf5_imagestacks import create_hdf5_imagestacks
from ..h5stacks.get_hdf5_imagestacks import get_hdf5_imagestacks
from ..common import timing
from ..io import nexus
from ..io import edf
from ..instruments import configuration

from .proc_math import execute as math
from .proc_align import execute as align
from .proc_replacevalue import execute as replacevalue
from .proc_crop import execute as execcrop
from .proc_common import defaultstack
from .proc_common import flattenstacks
from .proc_resample import execute as execresample

logger = logging.getLogger(__name__)

def getinstrument(kwargs):
    if "instrument" not in kwargs:
        raise RuntimeError("You need to specify an instrument.")
    instrument = kwargs["instrument"]
    if isinstance(instrument,configuration.InstrumentInfo):
        return instrument
    return configuration.factory(instrument,**kwargs.get("instrument_params",{}))

def exportedf(h5name,**kwargs):
    logger.info("EDF export {}:".format(h5name))
    
    instrument = getinstrument(kwargs)

    path = os.path.join(os.path.dirname(h5name),"results")
    if not os.path.exists(path):
        os.makedirs(path)

    filename = os.path.splitext(os.path.basename(h5name))[0]

    stacks, axes = get_hdf5_imagestacks(h5name,instrument.h5stackgroups)

    counters = instrument.counters()

    with h5py.File(h5name) as hdf5FileObject:
        for g in stacks:
            if "detectorsum" in stacks:
                if g.startswith("detector") and g!="detectorsum":
                    continue
        
            for s in stacks[g]:
                if s in counters:
                    continue
                energy = hdf5FileObject[g][s][instrument.edfheaderkeys["energylabel"]]
                n = len(energy)
                for i in range(n):
                    outfile = s.split("/")[-1]
                    outfile = outfile.replace("-","_")
                    outfile = outfile.replace("(","")
                    outfile = outfile.replace(")","")
                    if n==1:
                        outfile = os.path.join(path,"{}.edf".format(outfile))
                    else:
                        outfile = os.path.join(path,"{}_{}{}.edf".format(outfile,energy[i],instrument.edfheaderkeys["energyunit"]))

                    logger.info(outfile)
                    edf.saveedf(outfile,np.squeeze(hdf5FileObject[g][s]["data"][...,i]),{'Title': s},overwrite=True)

def createconfig_pre(sourcepath,destpath,scanname,scannumbers,cfgfiles,**kwargs):
    
    instrument = getinstrument(kwargs)
    stackdim = kwargs.get("stackdim",2)
    dtcor = kwargs.get("dtcor",True)
    mlines = kwargs.get("mlines",{})
    exclude_detectors = kwargs.get("exclude_detectors",None)
    include_detectors = kwargs.get("include_detectors",None)
    noxia = kwargs.get("noxia",False)
    addbeforefit = kwargs.get("addbeforefit",True)
    encodercor = kwargs.get("encodercor",{})
    fluxmonitor = kwargs.get("fluxmonitor",None)
    fluxmonitorparams = {"quantify":fluxmonitor is not None}
    counters = kwargs.get("counters",[])

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

    if not counters:
        lst = ["I0","It","If","calc"]
        if not noxia:
            lst.extend(["xrficr","xrfocr","xrfroi"])
        if encodercor:
            lst.extend(["motors"])
        counters = instrument.counters(include=lst)

    # Metadata is extracted from the counter headers
    # If they are missing, the xia stat files are checked
    metacounters = counters+["xia"]

    # Correct for deadtime when a single detector? (ignored when dtcor==False)
    # This exists because for one detector you can apply the deadtime correction
    # after XRF fitting
    dtcorifsingle = not all(k in instrument.counterdict for k in ["xrficr","xrfocr"])

    config = {
            # Input
            "sourcepath": sourcepath,
            "counter_reldir": instrument.counter_reldir,
            "scanname": scanname,
            "scannumbers": scannumbers,
            "counters": counters,
            "fluxcounter": instrument.counterdict["I0"],
            "transmissioncounter": instrument.counterdict["It"],

            # Meta data
            "metacounters": metacounters,
            "stacklabel": instrument.edfheaderkeys["energylabel"],
            "speccmdlabel": instrument.edfheaderkeys["speclabel"],
            "fastlabel": instrument.edfheaderkeys["fastlabel"],
            "slowlabel": instrument.edfheaderkeys["slowlabel"],
            "timelabel": instrument.edfheaderkeys["timelabel"],
            "stackinfo": instrument.imagemotors,

            # Deadtime correction
            "dtcor": dtcor,
            "dtcorifsingle":dtcorifsingle,

            # Configuration for fitting
            "detectorcfg": cfgfiles,
            "mlines": mlines,
            "fit": bfit,
            "fastfitting": True,
            "addbeforefitting": addbeforefit, # sum spectra
            "addafterfitting": True, # sum fit results and detector counters
            "exclude_detectors": exclude_detectors,
            "include_detectors": include_detectors,
            "fluxmonitor": fluxmonitorparams,
            
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

def process(sourcepath,destpath,scanname,scannumbers,cfgfiles,**kwargs):

    T0 = timing.taketimestamp()

    # Parse parameters
    instrument = getinstrument(kwargs)
    # ... stack
    bsamefile = False
    stackdim = kwargs.get("stackdim",2)
    default = kwargs.get("default",None)
    # ... xrf
    dtcor = kwargs.get("dtcor",True)
    fluxmonitor = kwargs.get("fluxmonitor",None)
    # ... normalization
    prealignnormcounter = kwargs.get("prealignnormcounter",None)
    postalignnormcounter = kwargs.get("postalignnormcounter",None)
    # ... align
    alignmethod = kwargs.get("alignmethod",None)
    alignreference = kwargs.get("alignreference",None)
    refimageindex = kwargs.get("refimageindex",None)
    roialign = kwargs.get("roialign",None)
    cropalign = False
    cropafter = kwargs.get("crop",False)
    replacenan = kwargs.get("replacenan",False)
    plot = kwargs.get("plot",False)
    # ... other
    encodercor = kwargs.get("encodercor",False)
    skippre = kwargs.get("skippre",False)
    
    # Image stacks (counters + result of XRF fitting)
    preprocessingexists = False
    if skippre:
        h5file = os.path.join(destpath,scanname[0]+".h5")
        preprocessingexists = os.path.isfile(h5file)

    if preprocessingexists:
        stacks, axes = get_hdf5_imagestacks(h5file,instrumentinfo.h5stackgroups)
    else:
        logger.info("Creating image stacks ...")
        jsonfile, h5file = createconfig_pre(sourcepath,destpath,scanname,scannumbers,\
                                            cfgfiles,**kwargs)
        stacks, axes = create_hdf5_imagestacks(jsonfile,fluxmonitor=fluxmonitor)

        #stacks2, axes2 = get_hdf5_imagestacks(h5file,["counters","detectorsum"])
        #assert(axes == axes2)
        #assert(stacks == stacks2)

    lastextension = ""

    if "detectorsum" in stacks or fluxmonitor is not None:
        dtcor = False # done on the raw data

    # Convert stack dictionary to stack list
    stacks = flattenstacks(stacks)

    # Default group
    defaultstack(h5file,stacks,default)

    # Groups that don't change and need to be copied
    copygroups = ["stackinfo"]

    # Normalization
    if dtcor or prealignnormcounter is not None:
        skip = instrument.counters(include=["xrfroi"])

        # Create normalization expression
        if dtcor:
            icr = instrument.counterdict["xrficr"]
            ocr = instrument.counterdict["xrfocr"]
            
            if prealignnormcounter is None:
                expression = "{{}}*nanone({{{}}}/{{{}}})".format(icr,ocr)
            else:
                expression = "{{}}*nanone({{{}}}/({{{}}}*{{{}}}))".format(icr,ocr,prealignnormcounter)
            skip += [icr,ocr]
        else:
            expression = "{{}}/{{{}}}".format(prealignnormcounter)

        if prealignnormcounter is not None:
            skip += [prealignnormcounter]

        h5file,stacks,axes = math(h5file,stacks,axes,copygroups,bsamefile,default,expression,skip,stackdim=stackdim,extension="norm")
        lastextension = "norm"
        
    # Correct for encoder positions
    if encodercor and instrument.encoderresolution:
        resampleinfo = {mot:{"encoder":instrument.counterdict["encoders"][mot],"resolution":res} for mot,res in instrument.encoderresolution.items()}
        if resampleinfo:
            h5file,stacks,axes = execresample(h5file, stacks, axes, copygroups, bsamefile, default, resampleinfo,extension="resample")
            lastextension = "resample"
            
    # Alignment
    if alignmethod is not None and alignreference is not None:
        # Alignment
        h5file,stacks,axes = align(h5file,stacks,axes, copygroups, bsamefile, default,\
            alignmethod, alignreference, refimageindex, cropalign, roialign, plot, stackdim)
        lastextension = "align"
        
    # Post normalization
    if postalignnormcounter is not None:
        skip = instrument.counters(include=["xrfroi"])
        skip.append(postalignnormcounter)
        skip.extend(["calc_flux0","calc_fluxt"])
        expression = "{{}}/{{{}}}".format(postalignnormcounter)
        h5file,stacks,axes = math(h5file,stacks,axes,copygroups,bsamefile,default,expression,skip,stackdim=stackdim,extension="postnorm")
        lastextension = "postnorm"
        
    # Remove NaN's
    if replacenan:
        orgvalue = np.nan
        newvalue = 0
        replacevalue(h5file, stacks, axes, copygroups, bsamefile, default, orgvalue, newvalue)
        lastextension = "replace"
        
    # Crop
    if cropafter:
        cropinfo = {"nanval":np.nan,"stackdim":stackdim,"reference set":alignreference}
        execcrop(h5file, stacks, axes, copygroups, bsamefile, default, cropinfo)
        #lastextension = "crop"
          
    timing.printtimeelapsed(T0,logger)
    
    return lastextension
    

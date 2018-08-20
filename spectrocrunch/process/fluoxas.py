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
import shutil

from ..xrf.create_hdf5_imagestacks import create_hdf5_imagestacks
from ..h5stacks.get_hdf5_imagestacks import get_hdf5_imagestacks
from ..common import timing
from ..io import nexus
from ..io import edf
from ..io.utils import mkdir
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
    return configuration.factory(instrument,**kwargs.get("instrument_parameters",{}))

def exportedf(h5name,**kwargs):
    logger.info("EDF export {}:".format(h5name))
    
    instrument = getinstrument(kwargs)

    path = os.path.basename(h5name)
    n = 0
    while len(path)!=n:
        n = len(path)
        path = os.path.splitext(path)[0]
    path = os.path.join(os.path.dirname(h5name),"{}_results".format(path))
    
    # not necessary but clean in case of reruns
    if os.path.isdir(path):
        shutil.rmtree(path)
            
    filename = os.path.splitext(os.path.basename(h5name))[0]

    stacks, axes, procinfo = get_hdf5_imagestacks(h5name,instrument.h5stackgroups)
    counters = instrument.counters()

    if "detectorsum" in stacks:
        stacks = {g:stacks[g] for g in stacks if not (g.startswith("detector") and g!="detectorsum")}
    ndet = sum(1 for k in stacks if k.startswith("detector"))

    with h5py.File(h5name) as hdf5FileObject:
        for g in stacks:
            if ndet==1:
                outpath = path
            else:
                outpath = os.path.join(path,g)
            
            if not os.path.exists(outpath):
                os.makedirs(outpath)
        
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
                        outfile = os.path.join(outpath,"{}.edf".format(outfile))
                    else:
                        outfile = os.path.join(outpath,"{}_{}{}.edf".format(outfile,energy[i],instrument.edfheaderkeys["energyunit"]))

                    logger.info(outfile)
                    edf.saveedf(outfile,np.squeeze(hdf5FileObject[g][s]["data"][...,i]),{'Title': s},overwrite=True)

def createconfig_pre(sourcepath,destpath,scanname,scannumbers,cfgfiles,**kwargs):
    
    instrument = getinstrument(kwargs)
    stackdim = kwargs.get("stackdim",2)
    dtcor = kwargs.get("dtcor",True)
    fastfitting = kwargs.get("fastfitting",True)
    addbeforefit = kwargs.get("addbeforefit",True)
    addafterfitting = kwargs.get("addafterfitting",True)
    mlines = kwargs.get("mlines",{})
    exclude_detectors = kwargs.get("exclude_detectors",None)
    include_detectors = kwargs.get("include_detectors",None)
    noxia = kwargs.get("noxia",False)
    encodercor = kwargs.get("encodercor",{})
    qxrfgeometry = kwargs.get("qxrfgeometry",None)
    qxrfgeometryparams = {"quantify":qxrfgeometry is not None}
    counters = kwargs.get("counters",[])
    fluxid = kwargs.get("fluxid","I0")
    transmissionid = kwargs.get("transmissionid","It")

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
        #lst = ["I0_counts","It_counts","If_counts","calc"]
        #if not noxia:
        #    lst.extend(["xrficr","xrfocr","xrfroi"])
        #if encodercor:
        #    lst.extend(["motors"])
        #counters = instrument.counters(include=lst)
        lst = []
        if noxia:
            lst.extend(["xrficr","xrfocr","xrfroi"])
        if not encodercor:
            lst.extend(["motors"])
        counters = instrument.counters(exclude=lst)

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
            "fluxcounter": instrument.counterdict[fluxid+"_counts"],
            "transmissioncounter": instrument.counterdict[transmissionid+"_counts"],

            # Meta data
            "metadata": instrument.metadata,
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
            "fastfitting": fastfitting,
            "addbeforefitting": addbeforefit, # sum spectra
            "addafterfitting": addafterfitting, # sum fit results and detector counters
            "exclude_detectors": exclude_detectors,
            "include_detectors": include_detectors,
            "qxrfgeometry": qxrfgeometryparams,
            
            # Output directories
            "outdatapath": os.path.join(destpath,scanname[0]+"_data"),
            "outfitpath": os.path.join(destpath,scanname[0]+"_fit"),
            "hdf5output": os.path.join(destpath,scanname[0]+".h5"),
            "stackdim": stackdim
    }

    # Create configuration file
    mkdir(destpath)
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
    qxrfgeometry = kwargs.get("qxrfgeometry",None)
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
        stacks, axes, procinfo = get_hdf5_imagestacks(h5file,instrument.h5stackgroups)
    else:
        jsonfile, h5file = createconfig_pre(sourcepath,destpath,scanname,scannumbers,\
                                            cfgfiles,**kwargs)
        stacks, axes, procinfo = create_hdf5_imagestacks(jsonfile,qxrfgeometry=qxrfgeometry)

        #stacks2, axes2, procinfo2 = get_hdf5_imagestacks(h5file,["counters","detectorsum"])
        #assert(axes == axes2)
        #assert(stacks == stacks2)

    h5filelast = h5file

    # Convert stack dictionary to stack list
    stacks = flattenstacks(stacks)

    # Default group
    defaultstack(h5file,stacks,default)

    # Groups that don't change and need to be copied
    copygroups = ["stackinfo"]

    # Normalization
    dtcor = procinfo[-1]["dtneeded"]
    if dtcor or prealignnormcounter is not None:
        skip = instrument.counters(include=["counters"])

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
        h5filelast = h5file
        
    # Correct for encoder positions
    if encodercor and instrument.encoderresolution:
        resampleinfo = {mot:{"encoder":instrument.counterdict["encoders"][mot],"resolution":res} for mot,res in instrument.encoderresolution.items()}
        if resampleinfo:
            h5file,stacks,axes = execresample(h5file, stacks, axes, copygroups, bsamefile, default, resampleinfo,extension="resample")
            h5filelast = h5file
            
    # Alignment
    if alignmethod is not None and alignreference is not None:
        # Alignment
        h5file,stacks,axes = align(h5file,stacks,axes, copygroups, bsamefile, default,\
            alignmethod, alignreference, refimageindex, cropalign, roialign, plot, stackdim)
        h5filelast = h5file
        
    # Post normalization
    if postalignnormcounter is not None:
        skip = instrument.counters(include=["xrfroi"])
        skip.append(postalignnormcounter)
        skip.extend(["calc_flux0","calc_fluxt"])
        expression = "{{}}/{{{}}}".format(postalignnormcounter)
        h5file,stacks,axes = math(h5file,stacks,axes,copygroups,bsamefile,default,expression,skip,stackdim=stackdim,extension="postnorm")
        h5filelast = h5file
        
    # Remove NaN's
    if replacenan:
        orgvalue = np.nan
        newvalue = 0
        h5filelast,_,_ = replacevalue(h5file, stacks, axes, copygroups, bsamefile, default, orgvalue, newvalue)
        
    # Crop
    if cropafter:
        cropinfo = {"nanval":np.nan,"stackdim":stackdim,"reference set":alignreference}
        h5filelast,_,_ = execcrop(h5file, stacks, axes, copygroups, bsamefile, default, cropinfo)
          
    timing.printtimeelapsed(T0,logger)
    
    return h5filelast
    

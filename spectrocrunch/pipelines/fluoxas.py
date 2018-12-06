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
import numpy as np
from ..utils import instance
from ..process import basetask
from ..io import nxfs
from ..instruments.configuration import getinstrument

def xrfparameters(**parameters):
    sourcepath = parameters["sourcepath"]
    scanname = parameters["scanname"]
    scannumbers = parameters["scannumbers"]
    cfgfiles = parameters["cfgfiles"]
    nxentry = parameters["nxentry"]
    
    instrument = getinstrument(parameters)
    dtcor = parameters.get("dtcor",True)
    fastfitting = parameters.get("fastfitting",True)
    adddetectors = parameters.get("adddetectors",True)
    addbeforefit = parameters.get("addbeforefit",True)
    mlines = parameters.get("mlines",{})
    exclude_detectors = parameters.get("exclude_detectors",None)
    include_detectors = parameters.get("include_detectors",None)
    noxia = parameters.get("noxia",False)
    encodercor = parameters.get("encodercor",{})
    qxrfgeometry = parameters.get("qxrfgeometry",None)
    correctspectra = parameters.get("correctspectra",False)
    fluxid = parameters.get("fluxid","I0")
    transmissionid = parameters.get("transmissionid","It")
    #dtcorcounters = all(k in instrument.counterdict for k in ["xrficr","xrfocr"])
    
    if noxia:
        cfgfiles = None
    bfit = cfgfiles is not None

    if exclude_detectors is None:
        exclude_detectors = []

    if include_detectors is None:
        include_detectors = []

    if not instance.isarray(sourcepath):
        sourcepath = [sourcepath]
    if not instance.isarray(scanname):
        scanname = [scanname]
    if not instance.isarray(scannumbers):
        scannumbers = [scannumbers]
    if not instance.isarray(cfgfiles):
        cfgfiles = [cfgfiles]

    lst = []
    if noxia:
        lst.extend(["xrficr","xrfocr","xrfroi"])
    if not encodercor:
        lst.extend(["motors"])
    counters = instrument.counters(exclude=lst)
    counters.extend(parameters.get("counters",[]))
    
    outputparent = nxfs.Path(str(nxentry))
    outputparent = outputparent.parent[outputparent.name+'.1']
    
    edffields1 = ('speclabel','slowlabel','fastlabel','energylabel','timelabel')
    edffields2 = ('speclabel','slowlabel','fastlabel','stackvalue','time')
    edfheader = {k2:instrument.edfheaderkeys[k1] for k1,k2 in zip(edffields1,edffields2)}
    edfheader['axesnamemap'] = instrument.imagemotors
    edfheader['compensationmotors'] = instrument.compensationmotors
    
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
            "edfheader": edfheader,
            "units": instrument.units,

            # Data correction
            "dtcor": dtcor,
            "correctspectra": correctspectra,
            "adddetectors": adddetectors, # sum spectra
            "addbeforefit": addbeforefit, # sum fit results and detector counters
            "qxrfgeometry": qxrfgeometry,
            
            # Configuration for fitting
            "detectorcfg": cfgfiles,
            "mlines": mlines,
            "fit": bfit,
            "fastfitting": fastfitting,
            "exclude_detectors": exclude_detectors,
            "include_detectors": include_detectors,

            # Output directories
            "outputparent": outputparent
    }
    
    return config,instrument
    
def tasks(**parameters):
    """Scripted pipeline to process FluoXAS data
    """
    tasks = []
    
    # Common parameters
    parameters['stackdim'] = parameters.get('stackdim',0)
    parameters['default'] = parameters.get('default',None)
    commonparams = {k:parameters[k] for k in ['default','stackdim']}

    # Image stacks (counters + result of XRF fitting)
    xrfparams,instrument = xrfparameters(**parameters)
    xrfparams.update(commonparams)
    task = basetask.task(method='pymca',name='process:pymca',**xrfparams)
    tasks.append(task)
    
    # Normalization
    prealignnormcounter = parameters.get("prealignnormcounter",None)
    dtcor = False # No longer needed
    if dtcor or prealignnormcounter is not None:
        copy = [{'method':'regexparent','pattern':prefix}
                for prefix in instrument.counterdict["counters"]]

        # Create normalization expression
        if dtcor:
            icr = instrument.counterdict["xrficr"]
            ocr = instrument.counterdict["xrfocr"]
            if prealignnormcounter is None:
                expression = "{{}}*nanone({{{}}}/{{{}}})".format(icr,ocr)
            else:
                expression = "{{}}*nanone({{{}}}/({{{}}}*{{{}}}))".format(icr,ocr,prealignnormcounter)
        else:
            expression = "{{}}/{{{}}}".format(prealignnormcounter)

        task = basetask.task(dependencies=task,method='expression',name='process:normalize',
                              expression=expression,copy=copy,**commonparams)
        tasks.append(task)
        
    # Correct for encoder positions
    encodercor = parameters.get("encodercor",False)
    if encodercor and instrument.encoderresolution:
        encoders = instrument.encoderinfo
        task = basetask.task(dependencies=task,method='resample',name='process:resample',
                              encoders=encoders,**commonparams)
        tasks.append(task)
            
    # Alignment
    alignmethod = parameters.get("alignmethod",None)
    alignreference = parameters.get("alignreference",None)
    if alignmethod and alignreference is not None:
        refimageindex = parameters.get("refimageindex",-1)
        roi = parameters.get("roialign",None)
        plot = parameters.get("plot",False)
        task = basetask.task(dependencies=task,method='align',name='process:align',alignmethod=alignmethod,
                              reference=alignreference,refimageindex=refimageindex,
                              crop=False,roi=roi,plot=plot,**commonparams)
        tasks.append(task)

    # Post normalization
    postalignnormcounter = parameters.get("postalignnormcounter",None)
    if postalignnormcounter is not None:
        copy = [{'method':'regexparent','pattern':prefix}
                for prefix in instrument.counterdict["counters"]]
        
        expression = "{{}}/{{{}}}".format(postalignnormcounter)
        task  = basetask.task(dependencies=task,method='expression',name='process:postnormalize',
                               expression=expression,copy=copy,**commonparams)
        tasks.append(task)
        
    # Remove NaN's
    replacenan = parameters.get("replacenan",False)
    if replacenan:
        tmp = basetask.task(dependencies=task,method='replace',name='process:replace',
                             org=np.nan,new=0,**commonparams)
        tasks.append(tmp)
                                            
    # Crop
    cropafter = parameters.get("crop",False)
    if cropafter:
        tmp = basetask.task(dependencies=task,method='crop',name='process:crop',nanval=np.nan,
                             reference=alignreference,**commonparams)
        tasks.append(tmp)
                                            
    return tasks

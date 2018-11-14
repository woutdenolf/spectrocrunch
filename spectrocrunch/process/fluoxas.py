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
import logging
import numpy as np

from ..xrf import create_hdf5_imagestacks
from ..utils import timing
from ..utils import instance
from ..data import nxtask
from ..io import nxfs
from ..instruments.configuration import getinstrument

logger = logging.getLogger(__name__)

def xrfparameters(**parameters):
    
    sourcepath = parameters["sourcepath"]
    destpath = parameters["destpath"]
    scanname = parameters["scanname"]
    scannumbers = parameters["scannumbers"]
    cfgfiles = parameters["cfgfiles"]
    
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
    nxentry = parameters.get("nxentry",None)
    
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
    
    h5file = os.path.join(destpath,scanname[0]+".h5")
    if nxentry:
        nxentry = nxfs.Path('/',h5file=h5file).nxentry(name=nxentry)
    else:
        nxentry = nxfs.Path('/',h5file=h5file).new_nxentry()
    
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
            "positionmotors": instrument.imagemotors,
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
            "destpath": destpath,
            "outbase": scanname[0],
            "outdatapath": os.path.join(destpath,scanname[0]+"_data"),
            "outfitpath": os.path.join(destpath,scanname[0]+"_fit"),
            "nxentry": nxentry
    }
    
    return config,instrument


def runtask(**parameters):
    task = nxtask.newtask(**parameters)
    task.run()
    nxprocess = task.output
    if nxprocess.exists:
        return nxprocess
    else:
        return None
        
def process(**parameters):

    with timing.timeit_logger(logger):
    
        # Common parameters
        parameters['stackdim'] = parameters.get('stackdim',0)
        parameters['default'] = parameters.get('default',None)
        commonparams = {k:parameters[k] for k in ['default','stackdim']}

        # Image stacks (counters + result of XRF fitting)
        xrfparams,instrument = xrfparameters(**parameters)
        xrfparams.update(commonparams)
        nxprocess = runtask(method='xrf',**xrfparams)
        if nxprocess is None:
            return nxprocess
        else:
            nxprocesslast = nxprocess
        
        # Normalization
        prealignnormcounter = parameters.get("prealignnormcounter",None)
        dtcor = not nxprocess.results['info']['dtneeded'] # No longer needed
        if dtcor or prealignnormcounter is not None:
            skip = [{'method':'regexparent','pattern':'counters'},
                    {'method':'regex','pattern':instrument.counterdict["xrficr"]},
                    {'method':'regex','pattern':instrument.counterdict["xrfocr"]},
                    ]

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

            nxprocess = runtask(previous=nxprocess,method='expression',name='normalize',
                                expression=expression,skip=skip,**commonparams)
            if nxprocess is None:
                return nxprocesslast
            else:
                nxprocesslast = nxprocess
            
        # Correct for encoder positions
        encodercor = parameters.get("encodercor",False)
        if encodercor and instrument.encoderresolution:
            encoders = instrument.encoderinfo
            nxprocess = runtask(previous=nxprocess,method='resample',
                                encoders=encoders,**commonparams)
            if nxprocess is None:
                return nxprocesslast
            else:
                nxprocesslast = nxprocess
                
        # Alignment
        alignmethod = parameters.get("alignmethod",None)
        alignreference = parameters.get("alignreference",None)
        if alignmethod and alignreference is not None:
            refimageindex = parameters.get("refimageindex",-1)
            roi = parameters.get("roialign",None)
            plot = parameters.get("plot",False)
            nxprocess = runtask(previous=nxprocess,method='align',alignmethod=alignmethod,
                                reference=alignreference,refimageindex=refimageindex,
                                crop=False,roi=roi,plot=plot,**commonparams)
            if nxprocess is None:
                return nxprocesslast
            else:
                nxprocesslast = nxprocess

        # Post normalization
        postalignnormcounter = parameters.get("postalignnormcounter",None)
        if postalignnormcounter is not None:
            skip = [{'method':'regexparent','pattern':'counters'},
                    {'method':'regex','pattern':instrument.counterdict["xrficr"]},
                    {'method':'regex','pattern':instrument.counterdict["xrfocr"]},
                    ]
            
            expression = "{{}}/{{{}}}".format(postalignnormcounter)
            nxprocess = runtask(previous=nxprocess,method='expression',name='postnormalize',
                                expression=expression,skip=skip,**commonparams)
            if nxprocess is None:
                return nxprocesslast
            else:
                nxprocesslast = nxprocess
            
        # Remove NaN's
        replacenan = parameters.get("replacenan",False)
        if replacenan:
            tmp = runtask(nxprocess,method='replace',
                          org=np.nan,new=0,**commonparams)
            if tmp is not None:
                nxprocesslast = tmp
                                                
        # Crop
        cropafter = parameters.get("crop",False)
        if cropafter:
            tmp = runtask(previous=nxprocess,method='crop',nanval=np.nan,
                          reference=alignreference,**commonparams)
            if tmp is not None:
                nxprocesslast = tmp
                                                
        return nxprocesslast

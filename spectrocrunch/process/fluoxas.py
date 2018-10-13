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
from ..h5stacks.get_hdf5_imagestacks import get_hdf5_imagestacks
from ..utils import timing
from ..utils import instance
from ..io import edf

from .proc_math import execute as math
from .proc_align import execute as align
from .proc_replacevalue import execute as replacevalue
from .proc_crop import execute as execcrop
from .proc_utils import defaultstack
from .proc_utils import flattenstacks
from .proc_resample import execute as execresample

logger = logging.getLogger(__name__)
        
def process(**parameters):

    T0 = timing.taketimestamp()

    # Image stacks (counters + result of XRF fitting)
    instrument = create_hdf5_imagestacks.getinstrument(parameters)
    default = parameters.get("default",None)
    skippre = parameters.get("skippre",False)
 
    h5file = create_hdf5_imagestacks.h5file(**parameters)
    if skippre:
        preprocessingexists = os.path.isfile(h5file)
    else:
        preprocessingexists = False
        
    if preprocessingexists:
        stacks, axes, procinfo = get_hdf5_imagestacks(h5file,instrument.h5stackgroups)
    else:
        stacks, axes, procinfo = create_hdf5_imagestacks.create_hdf5_imagestacks(**parameters)
    defaultstack(h5file,stacks,default)
    stacks = flattenstacks(stacks)
    h5filelast = h5file
    
    # Common parameters for the following steps
    bsamefile = False
    stackdim = parameters.get("stackdim",2)
    copygroups = ["stackinfo"]

    # Normalization
    prealignnormcounter = parameters.get("prealignnormcounter",None)
    dtcor = procinfo[-1]["result"]["dtneeded"] # No longer needed
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
    encodercor = parameters.get("encodercor",{})
    if encodercor and instrument.encoderresolution:
        resampleinfo = {mot:{"encoder":instrument.counterdict["encoders"][mot],"resolution":res} for mot,res in instrument.encoderresolution.items()}
        if resampleinfo:
            h5file,stacks,axes = execresample(h5file, stacks, axes, copygroups, bsamefile, default, resampleinfo,extension="resample")
            h5filelast = h5file
            
    # Alignment
    alignmethod = parameters.get("alignmethod",None)
    alignreference = parameters.get("alignreference",None)
    if alignmethod and alignreference is not None:
        refimageindex = parameters.get("refimageindex",None)
        roialign = parameters.get("roialign",None)
        cropalign = False
        plot = parameters.get("plot",False)
    
        # Alignment
        h5file,stacks,axes = align(h5file,stacks,axes, copygroups, bsamefile, default,\
            alignmethod, alignreference, refimageindex, cropalign, roialign, plot, stackdim)
        h5filelast = h5file

    # Post normalization
    postalignnormcounter = parameters.get("postalignnormcounter",None)
    if postalignnormcounter is not None:
        skip = instrument.counters(include=["xrfroi"])
        skip.append(postalignnormcounter)
        skip.extend(["calc_flux0","calc_fluxt"])
        expression = "{{}}/{{{}}}".format(postalignnormcounter)
        h5file,stacks,axes = math(h5file,stacks,axes,copygroups,bsamefile,default,expression,skip,stackdim=stackdim,extension="postnorm")
        h5filelast = h5file
        
    # Remove NaN's
    replacenan = parameters.get("replacenan",False)
    if replacenan:
        orgvalue = np.nan
        newvalue = 0
        h5filelast,_,_ = replacevalue(h5file, stacks, axes, copygroups, bsamefile, default, orgvalue, newvalue)
        
    # Crop
    cropafter = parameters.get("crop",False)
    if cropafter:
        cropinfo = {"nanval":np.nan,"stackdim":stackdim,"reference set":alignreference}
        h5filelast,_,_ = execcrop(h5file, stacks, axes, copygroups, bsamefile, default, cropinfo)
          
    timing.printtimeelapsed(T0,logger)
    
    return h5filelast
    

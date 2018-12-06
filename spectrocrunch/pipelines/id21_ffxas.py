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

import os
from ..utils import instance
from ..process import basetask
from ..io import nxfs

def ffparameters(**parameters):
    sourcepath = parameters["sourcepath"]
    radix = parameters["radix"]

    rebin = parameters.get("rebin",(1,1))
    roiraw = parameters.get("roiraw",None)
    stackdim = parameters.get("stackdim",2)
    normalize = parameters.get("normalize",False)
    normalizeonload = parameters.get("normalizeonload",True)
    flatbeforeafter = parameters.get("flatbeforeafter",True)
    nxentry = parameters.get("nxentry",None)
    
    if not instance.isarray(sourcepath):
        sourcepath = [sourcepath]
    if not instance.isarray(radix):
        radix = [radix]

    outputparent = nxfs.Path(str(nxentry))
    outputparent = outputparent.parent[outputparent.name+'.1']

    config = {
        # EDF header
        "frametimelabel" : "exposure_time",
        "frametimedefault" : 0.,
        "roilabel" : "roi",
        "nbframeslabel" : "nb_frames",
        "stacklabel" : "energy",

        # Defaults
        "dtype": "np.float32", # must be floating point!
        "darkcurrentzero" : 90.,
        "darkcurrentgain" : 1.,
        "stackdim" : stackdim,

        # Data
        "darklist" : map(lambda xy: os.path.join(xy[0],xy[1]+"*_dark_*.edf"),zip(sourcepath,radix)),
        "datalist" : map(lambda xy: os.path.join(xy[0],xy[1]+"*_data_*.edf"),zip(sourcepath,radix)),
        "flatlist" : map(lambda xy: os.path.join(xy[0],xy[1]+"*_ref_*.edf"),zip(sourcepath,radix)),
        "flatbeforeafter" : flatbeforeafter, # split up flat fields in before and after (same number of images)
        "normalize": normalize and normalizeonload,
        "keepflat": normalize and not normalizeonload,
        "roi": roiraw,
        "rebin":rebin,

        # Output
        "outputparent": outputparent,
        }
        
    return config

def tasks(**parameters):
    tasks = []
    
    # Common parameters
    parameters['stackdim'] = parameters.get('stackdim',0)
    parameters['default'] = parameters.get('default','sample')
    commonparams = {k:parameters[k] for k in ['default','stackdim']}

    # Image stacks (fullfield)
    ffparams = ffparameters(**parameters)
    ffparams.update(commonparams)
    task = basetask.task(method='fullfield',name='process:fullfield',**ffparams)
    tasks.append(task)
    
    # Normalization
    normalize = parameters.get("normalize",False)
    if normalize and not ffparams['normalize']:
        skip = [{'method':'regex','pattern':'flat1'},
                {'method':'regex','pattern':'flat2'},
                ]
        if ffparams['flatbeforeafter']:
            expression = "-ln(2*{}/({flat1}+{flat2}))"
        else:
            expression = "-ln({}/{flat1})"
        task = basetask.task(dependencies=task,method='expression',name='process:normalize',
                              expression=expression,skip=skip,**commonparams)
        tasks.append(task)
    
    # Alignment
    alignmethod = parameters.get("alignmethod",None)
    alignreference = 'sample'
    if alignmethod and alignreference is not None:
        refimageindex = parameters.get("refimageindex",-1)
        roi = parameters.get("roialign",None)
        plot = parameters.get("plot",False)
        task = basetask.task(dependencies=task,method='align',name='process:align',alignmethod=alignmethod,
                              reference=alignreference,refimageindex=refimageindex,
                              crop=False,roi=roi,plot=plot,**commonparams)
        tasks.append(task)
        
    # Crop
    roiresult = parameters.get("roiresult",None)
    if roiresult:
        tmp = basetask.task(dependencies=task,method='crop',name='process:roi',roi=roiresult,
                             reference=alignreference,**commonparams)
        tasks.append(tmp)

    return tasks

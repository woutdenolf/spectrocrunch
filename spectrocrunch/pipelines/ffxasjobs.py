# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
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
from . import batch
from .run import run_sequential
from .id21_ffxas import tasks as ffxas_tasks
from ..utils import instance
from ..instruments.configuration import getinstrument

def staticscan(samplename,datasetname,radix,**parameters):
    jobname = batch.jobname("static",(samplename,datasetname,radix),parameters)

    instrument = getinstrument(parameters)
    mradix,subdir = instrument.fflocation(samplename,datasetname,type="static")
    if not instance.isarray(radix):
        radix = [radix]
    sourcepath = [os.path.join(parameters["proposaldir"],subdir,rdx) for rdx in radix]
    
    if len(radix)>1:
        nxentry = '{}_{}'.format(mradix,radix[0],radix[-1])
    else:
        nxentry = '{}_{}'.format(mradix,radix[0])
    nxentry = os.path.join(parameters.get('resultsdir',''),mradix,mradix+'.h5')+':/'+nxentry
    processdata(jobname,sourcepath,radix,nxentry,**parameters)
    
def processdata(jobname,*args,**kwargs):
    if "jobs" in kwargs:
        kwargs["jobs"].append((jobname,processdata_exec,args,kwargs))
    else:
        processdata_exec(*args,**kwargs)
        
def processdata_exec(sourcepath,radix,nxentry,**kwargs):
    parameters = dict(kwargs)
    parameters['sourcepath'] = sourcepath
    parameters['radix'] = radix
    parameters['nxentry'] = nxentry
    tasks = ffxas_tasks(**parameters)
    if run_sequential(tasks,name='fullfield'):
        pass
    else:
        unfinished = [task for task in tasks if not task.done]
        raise RuntimeError('The following tasks are not finished: {}'.format(unfinished))

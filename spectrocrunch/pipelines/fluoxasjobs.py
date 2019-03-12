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
import logging

from . import batch
from .run import run_sequential
from . import fluoxas as mfluoxas
from ..process import nxresult
from ..utils import instance
from ..instruments.configuration import getinstrument
from ..io.edf import saveedf

logger = logging.getLogger(__name__)


def ensure_geometry_outputparent(parameters):
    params = mfluoxas.task_parameters(parameters, 'geometry')
    if params.get('outputparent', None):
        return
    nxroot = os.path.join(parameters.get('resultsdir', ''),
                          'geometries.h5')
    params['outputparent'] = nxroot+'::/xrf'


def ensure_outputparent(parameters, filename, entryname, scannumbers):
    ensure_geometry_outputparent(parameters)
    params = mfluoxas.task_parameters(parameters, 'common')
    if params.get('outputparent', None):
        return
    nxroot = os.path.join(parameters.get('resultsdir', ''),
                          filename+'.h5')
    if len(scannumbers) > 1:
        a, b = scannumbers[0], scannumbers[-1]
    else:
        a, b = scannumbers[0], None
    if instance.isarray(a):
        a = a[0]
    if instance.isarray(b):
        b = b[-1]
    if b:
        nxentry = "{}{}_{}".format(entryname, a, b)
    else:
        nxentry = "{}{}".format(entryname, a)
    params['outputparent'] = nxroot+'::/'+nxentry


def fluoxas(samplename, datasetname, scannumbers, mapnumbers, cfgfiles,
            **parameters):
    if len(scannumbers) != len(mapnumbers):
        raise RuntimeError(
            "fluoXAS map numbers must be equal to the mapnumbers")
    jobname = batch.jobname(
        "fluoxas",
        (samplename, datasetname, scannumbers, mapnumbers, cfgfiles),
        parameters)

    params = mfluoxas.task_parameters(parameters, 'common')
    instrument = getinstrument(**params)
    radix, subdir = instrument.xrflocation(
        samplename, datasetname, type="dynamic")

    sourcepaths = [os.path.join(parameters["proposaldir"], subdir, 
                   "{}_fluoXAS_{}".format(radix, nr)) for nr in scannumbers]
    scannames = ["{}_fluoXAS_{}".format(radix, nr) for nr in scannumbers]

    ensure_outputparent(parameters, samplename, radix+'.fluoxas', scannumbers)
    processdata(jobname, sourcepaths, scannames, mapnumbers,
                cfgfiles, fluoxas=True, **parameters)


def multi(samplename, datasetname, mapnumbers, cfgfiles, **parameters):
    jobname = batch.jobname(
        "multi",
        (samplename, datasetname, mapnumbers, cfgfiles),
        parameters)

    params = mfluoxas.task_parameters(parameters, 'common')
    instrument = getinstrument(**params)
    radix, subdir = instrument.xrflocation(
        samplename, datasetname, type="dynamic")
    sourcepaths = [os.path.join(parameters["proposaldir"], subdir)]

    scannames = [radix]
    scannumbers = [mapnumbers]
    ensure_outputparent(parameters, samplename, radix+'.sixes', mapnumbers)
    processdata(jobname, sourcepaths, scannames, scannumbers,
                cfgfiles, multi=True, **parameters)


def single(samplename, datasetname, mapnumber, cfgfiles, **parameters):
    jobname = batch.jobname(
        "single",
        (samplename, datasetname, mapnumber, cfgfiles),
        parameters)

    params = mfluoxas.task_parameters(parameters, 'common')
    instrument = getinstrument(**params)
    radix, subdir = instrument.xrflocation(
        samplename, datasetname, type="dynamic")
    sourcepaths = [os.path.join(parameters["proposaldir"], subdir)]

    scannames = [radix]
    scannumbers = [[mapnumber]]
    ensure_outputparent(parameters, samplename, radix+'.map', [mapnumber])
    processdata(jobname, sourcepaths, scannames,
                scannumbers, cfgfiles, **parameters)


def manualselection(sourcepaths, scannames, scannumbers, cfgfiles,
                    outname=None, outsuffix=None, **parameters):
    jobname = batch.jobname(
        "manualselection",
        (sourcepaths, scannames, scannumbers, cfgfiles),
        parameters)

    if outname is None:
        outname = scannames[0]
    if outsuffix is None:
        outsuffix = '.map'
    ensure_outputparent(parameters, outname, outname+outsuffix, scannumbers)
    processdata(jobname, sourcepaths, scannames,
                scannumbers, cfgfiles, **parameters)


def processdata(jobname, *args, **kwargs):
    jobs = kwargs.pop('jobs', None)
    if jobs is None:
        # Execute immediately
        processdata_exec(*args, **kwargs)
    else:
        # Add to queue
        jobs.append((jobname, processdata_exec, args, kwargs))


def processdata_exec(sourcepaths, scannames, scannumbers, cfgfiles,
                     fluoxas=False, multi=False, resultsdir=None,
                     edfexport=False, **parameters):
    # Basic input
    params = mfluoxas.task_parameters(parameters, 'pymca')
    params['sourcepaths'] = sourcepaths
    params['scannames'] = scannames
    params['scannumbers'] = scannumbers
    if not resultsdir:
        resultsdir = ''
    if instance.isstring(cfgfiles):
        cfgfiles = [cfgfiles]
    params['pymcacfg'] = [cfg if os.path.isabs(cfg)
                          else os.path.join(resultsdir, cfg)
                          for cfg in cfgfiles]

    # Image aligment
    params = mfluoxas.task_parameters(parameters, 'align')
    if not fluoxas and not multi:
        params['alignmethod'] = None

    # Process
    tasks = mfluoxas.tasks(**parameters)
    if edfexport:
        edfoutput = not fluoxas and not tasks[-1].done
    else:
        edfoutput = False
    if run_sequential(tasks, name='fluoxas'):
        if edfoutput:
            exportedf(tasks[-1].output)
    else:
        unfinished = [task for task in tasks if not task.done]
        raise RuntimeError(
            'The following tasks are not finished: {}'.format(unfinished))


def exportedf(nxprocess):
    outdir = nxprocess.device.parent[nxprocess.device.name+'_edfresults']
    logger.info("EDF export:\n Input: {}\n Output: {}".format(
        nxprocess, outdir))
    outdir.remove(recursive=True)
    outdir.mkdir()

    groups, axes, stackdim = nxresult.regulargriddata(nxprocess)
    stackaxes = axes[stackdim]
    for group, paths in groups.items():
        if group.isdetector:
            for path in paths:
                with path.open(mode='r') as dset:
                    shape = dset.shape
                    index = [slice(None)]*dset.ndim
                    if shape[stackdim] == 1:
                        filename = group.xialabel+'_'+path.name+'.edf'
                    else:
                        filename = group.xialabel+'_'+path.name+'_{}keV.edf'
                    for i in range(shape[stackdim]):
                        index[stackdim] = i
                        image = dset[tuple(index)]
                        name = filename.format(stackaxes[i])
                        title = '{}@{}keV'.format(group.xialabel, stackaxes[i])
                        logger.info(' saving {}'.format(filename))
                        saveedf(outdir[name].path, image, {
                                'Title': title}, overwrite=True)

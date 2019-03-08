# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the 'Software'), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import os
import numpy as np
from ..utils import instance
from ..process.utils import create_task
from ..io import nxfs
from ..instruments.configuration import getinstrument


def task_parameters(parameters, name):
    params = parameters.get(name, {})
    parameters[name] = params
    return params


def ensure_parameter(parameters, name, default=None):
    value = parameters.get(name, default)
    parameters[name] = value
    return value


def require_parameter(parameters, name):
    if name not in parameters:
        raise RuntimeError('Parameter {} is missing'.format(repr(name)))
    return parameters[name]


def require_pop(parameters, name):
    require_parameter(parameters, name)
    return parameters.pop(name)


def ensure_list(parameters, name, default=None):
    value = parameters.get(name, default)
    if value is None:
        value = []
    elif not instance.isarray(value):
        value = [value]
    parameters[name] = value
    return value


def require_list(parameters, name):
    require_parameter(parameters, name)
    return ensure_list(parameters, name)


def require_poplist(parameters, name):
    require_parameter(parameters, name)
    lst = ensure_list(parameters, name)
    parameters.pop(name)
    return lst


def xrfparameters(parameters, instrument, include_encoders=True, quant=False):
    # Input
    require_list(parameters, 'sourcepaths')
    require_list(parameters, 'scannames')
    require_list(parameters, 'scannumbers')
    require_list(parameters, 'pymcacfg')

    # Extract relevant instrument info
    excounters = []
    if include_encoders:
        excounters.extend(['motors'])
    #dtcorcounters = all(k in instrument.counterdict for k in ['xrficr','xrfocr'])
    #if nospectra:
    #    excounters.extend(['xrficr', 'xrfocr', 'xrfroi'])
    counters = instrument.counters(exclude=excounters)
    counters.extend(parameters.get('counters', []))
    parameters['counters'] = counters

    edffields1 = ('speclabel', 'slowlabel', 'fastlabel',
                  'energylabel', 'timelabel')
    edffields2 = ('speclabel', 'slowlabel', 'fastlabel', 'stackvalue', 'time')
    edfheader = {k2: instrument.edfheaderkeys[k1]
                 for k1, k2 in zip(edffields1, edffields2)}
    edfheader['axesnamemap'] = instrument.imagemotors
    edfheader['compensationmotors'] = instrument.compensationmotors
    parameters['edfheader'] = edfheader

    parameters['counter_reldir'] = instrument.counter_reldir
    fluxid = parameters.pop('fluxid', 'I0')
    transmissionid = parameters.pop('transmissionid', 'It')
    parameters['fluxcounter'] = instrument.counterdict[fluxid+'_counts']
    parameters['transmissioncounter'] = instrument.counterdict[transmissionid+'_counts']
    parameters['metadata'] = instrument.metadata
    parameters['units'] = instrument.units

    # Others
    ensure_parameter(parameters, 'dtcor', True)
    ensure_parameter(parameters, 'fastfitting', True)
    ensure_parameter(parameters, 'adddetectors', True)  # sum spectra
    ensure_parameter(parameters, 'addbeforefit', True)  # sum fit results and detector counters
    ensure_parameter(parameters, 'mlines', {})
    ensure_parameter(parameters, 'correctspectra', False)
    ensure_list(parameters, 'exclude_detectors')
    ensure_list(parameters, 'include_detectors')
    ensure_list(parameters, 'samplecovers')
    ensure_list(parameters, 'transmissionfilters')
    if quant:
        require_parameter(parameters, 'diodeI0gain')
        require_parameter(parameters, 'diodeItgain')
        require_list(parameters, 'xrf_positions')


def tasks(**parameters):
    """Scripted pipeline to process FluoXAS data
    """
    tasks = []

    # Common parameters
    commonparams = task_parameters(parameters, 'common').copy()
    instrument = require_pop(commonparams, 'instrument')
    outputparent = nxfs.Path(str(require_pop(commonparams, 'nxentry')))
    ensure_parameter(commonparams, 'stackdim', 0)
    ensure_parameter(commonparams, 'default', None)
    instrument = getinstrument(instrument=instrument)

    # Calibrate geometry
    params = task_parameters(parameters, 'geometry')
    geometry = ensure_parameter(params, 'geometry', None)
    if geometry:
        init = ensure_parameter(params, 'init', {})
        init['instrument'] = instrument
        task = create_task(method='xrfgeometry',
                           name='xrfgeometry',
                           outputparent=outputparent,
                           **params)
        outputparent = None
        tasks.append(task)
    else:
        task = None

    # XRF fit and counter maps
    params = task_parameters(parameters, 'resample')
    encodercor = ensure_parameter(params, 'encodercor', False)
    params = task_parameters(parameters, 'pymca').copy()
    xrfparameters(params, instrument, include_encoders=encodercor, quant=bool(geometry))
    params.update(commonparams)
    task = create_task(method='pymca', name='pymca',
                       dependencies=task, 
                       outputparent=outputparent, **params)
    tasks.append(task)

    # Normalization
    params = task_parameters(parameters, 'prealignnormalize')
    counter = ensure_parameter(params, 'counter', None)
    dtcor = False  # No longer needed
    if dtcor or counter:
        copy = [{'method': 'regex', 'pattern': prefix}
                for prefix in instrument.counterdict['counters']]
        # Create normalization expression
        if dtcor:
            icr = instrument.counterdict['xrficr']
            ocr = instrument.counterdict['xrfocr']
            if counter:
                expression = '{{}}*nanone({{{}}}/({{{}}}*{{{}}}))'.format(
                    icr, ocr, counter)
            else:
                expression = '{{}}*nanone({{{}}}/{{{}}})'.format(icr, ocr)
        else:
            expression = '{{}}/{{{}}}'.format(counter)
        task = create_task(dependencies=task,
                           method='expression',
                           name='normalize',
                           expression=expression,
                           copy=copy,
                           **commonparams)
        tasks.append(task)

    # Correct for encoder positions
    if encodercor:
        encoders = instrument.encoderinfo
        task = create_task(dependencies=task,
                           method='resample',
                           name='resample',
                           encoders=encoders,
                           **commonparams)
        tasks.append(task)

    # Alignment
    params = task_parameters(parameters, 'align')
    alignmethod = ensure_parameter(params, 'alignmethod', None)
    reference = ensure_parameter(params, 'reference', None)
    if alignmethod and reference is not None:
        ensure_parameter(params, 'refimageindex', -1)
        ensure_parameter(params, 'roi', None)
        ensure_parameter(params, 'plot', False)
        ensure_parameter(params, 'crop', False)
        params.update(commonparams)
        task = create_task(dependencies=task,
                           method='align',
                           name='align',
                           alignmethod=alignmethod,
                           reference=reference,
                           **params)
        tasks.append(task)

    # Post normalization
    params = task_parameters(parameters, 'prostalignnormalize')
    counter = ensure_parameter(params, 'counter', None)
    if counter:
        copy = [{'method': 'regexparent', 'pattern': prefix}
                for prefix in instrument.counterdict['counters']]
        expression = '{{}}/{{{}}}'.format(counter)
        task = create_task(dependencies=task, method='expression',
                           name='postnormalize', expression=expression,
                           copy=copy, **commonparams)
        tasks.append(task)

    # Remove NaN's
    params = task_parameters(parameters, 'replacenan')
    replacenan = ensure_parameter(params, 'replacenan', False)
    if replacenan:
        tmp = create_task(dependencies=task, method='replace',
                          name='replace', org=np.nan, new=0,
                          **commonparams)
        tasks.append(tmp)

    # Crop
    params = task_parameters(parameters, 'crop')
    cropafter = params.get('crop', False)
    if cropafter:
        cropfull = params.get('cropfull', True)
        tmp = create_task(dependencies=task, method='crop',
                          name='crop', nanval=np.nan,
                          nanfull=cropfull, reference=reference,
                          **commonparams)
        tasks.append(tmp)

    return tasks

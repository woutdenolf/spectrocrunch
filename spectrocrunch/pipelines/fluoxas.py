# -*- coding: utf-8 -*-

import os
import numpy as np
from ..utils import instance
from ..process.utils import create_task
from ..io import nxfs
from ..instruments.configuration import getinstrument


def task_parameters(parameters, name):
    params = parameters.get(name, {}).copy()
    parameters[name] = params
    return params


def ensure_parameter(parameters, name, default=None):
    value = parameters.get(name, default)
    parameters[name] = value
    return value


def require_parameter(parameters, name):
    if name not in parameters:
        raise RuntimeError("Parameter {} is missing".format(repr(name)))
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
    require_list(parameters, "sourcepaths")
    require_list(parameters, "scannames")
    require_list(parameters, "scannumbers")
    require_list(parameters, "pymcacfg")

    # Extract relevant instrument info
    excounters = []
    if include_encoders:
        excounters.extend(["motors"])
    # dtcorcounters = all(k in instrument.counterdict for k in ['xrficr','xrfocr'])
    # if nospectra:
    #    excounters.extend(['xrficr', 'xrfocr', 'xrfroi'])
    counters = instrument.counters(exclude=excounters)
    counters.extend(parameters.get("counters", []))
    parameters["counters"] = counters

    edffields1 = ("speclabel", "slowlabel", "fastlabel", "energylabel", "timelabel")
    edffields2 = ("speclabel", "slowlabel", "fastlabel", "stackvalue", "time")
    edfheader = {
        k2: instrument.edfheaderkeys[k1] for k1, k2 in zip(edffields1, edffields2)
    }
    edfheader["axesnamemap"] = instrument.imagemotors
    edfheader["compensationmotors"] = instrument.compensationmotors
    parameters["edfheader"] = edfheader

    parameters["counter_reldir"] = instrument.counter_reldir
    fluxid = parameters.pop("fluxid", "I0")
    transmissionid = parameters.pop("transmissionid", "It")
    parameters["fluxcounter"] = instrument.counterdict[fluxid + "_counts"]
    parameters["transmissioncounter"] = instrument.counterdict[
        transmissionid + "_counts"
    ]
    parameters["metadata"] = instrument.metadata
    parameters["units"] = instrument.units

    # Others
    ensure_parameter(parameters, "dtcor", True)
    ensure_parameter(parameters, "fastfitting", True)
    ensure_parameter(parameters, "adddetectors", True)  # sum spectra
    ensure_parameter(
        parameters, "addbeforefit", True
    )  # sum fit results and detector counters
    ensure_parameter(parameters, "mlines", {})
    ensure_parameter(parameters, "addhigh", 0)
    ensure_parameter(parameters, "correctspectra", False)
    ensure_list(parameters, "exclude_detectors")
    ensure_list(parameters, "include_detectors")
    ensure_list(parameters, "samplecovers")
    ensure_list(parameters, "transmissionfilters")
    if quant and False:
        require_parameter(parameters, "gaindiodeI0")
        require_parameter(parameters, "gaindiodeIt")
        require_list(parameters, "xrf_positions")


def tasks(**parameters):
    """Scripted pipeline to process FluoXAS data"""
    tasks = []

    # Common parameters
    commonparams = task_parameters(parameters, "common")
    sinstrument = require_pop(commonparams, "instrument")
    coutputparent = commonparams.pop("outputparent", None)
    ensure_parameter(commonparams, "stackdim", 0)
    ensure_parameter(commonparams, "default", None)
    instrument = getinstrument(instrument=sinstrument)

    # Calibrate geometry
    params = task_parameters(parameters, "geometry")
    ensure_parameter(params, "outputparent", default=None)
    geometry = ensure_parameter(params, "geometry", None)
    if geometry:
        init = ensure_parameter(params, "init", {})
        init["instrument"] = instrument
        task = create_task(method="xrfgeometry", name="xrfgeometry", **params)
        tasks.append(task)
    else:
        task = None

    # XRF fit and counter maps
    params = task_parameters(parameters, "resample")
    encodercor = params.pop("encodercor", False)
    params = task_parameters(parameters, "pymca")
    ensure_parameter(params, "outputparent", default=coutputparent)
    xrfparameters(params, instrument, include_encoders=encodercor, quant=bool(geometry))
    params.update(commonparams)
    task = create_task(method="pymca", name="pymca", dependencies=task, **params)
    tasks.append(task)

    # Normalization
    params = task_parameters(parameters, "prealignnormalize")
    counter = params.pop("counter", None)
    if counter:
        copy = [
            {"method": "regex", "pattern": prefix}
            for prefix in instrument.counterdict["counters"]
        ]
        # Create normalization expression
        target = params.pop("target", None)
        if not target:
            target = f"mean({{{counter}}})"
        expression = f"{{}}*{target}/{{{counter}}}"
        ensure_parameter(params, "outputparent", default=coutputparent)
        ensure_parameter(params, "name", "prenormalize")
        task = create_task(
            dependencies=task,
            method="expression",
            name="normalize",
            expression=expression,
            copy=copy,
            **commonparams,
        )
        tasks.append(task)

    # Correct for encoder positions
    if encodercor:
        params = task_parameters(parameters, "resample")
        params.update(commonparams)
        ensure_parameter(params, "outputparent", default=coutputparent)
        ensure_parameter(params, "name", "resample")
        encoders = instrument.encoderinfo
        task = create_task(
            dependencies=task, method="resample", encoders=encoders, **params
        )
        tasks.append(task)

    # Alignment
    params = task_parameters(parameters, "align")
    alignmethod = ensure_parameter(params, "alignmethod", None)
    reference = ensure_parameter(params, "reference", None)
    if alignmethod and reference is not None:
        params.update(commonparams)
        ensure_parameter(params, "outputparent", default=coutputparent)
        ensure_parameter(params, "name", "align")
        ensure_parameter(params, "refimageindex", -1)
        ensure_parameter(params, "roi", None)
        ensure_parameter(params, "plot", False)
        ensure_parameter(params, "crop", False)
        ensure_parameter(params, "name", "postnormalize")
        task = create_task(dependencies=task, method="align", **params)
        tasks.append(task)

    # Post normalization
    params = task_parameters(parameters, "postalignnormalize")
    counter = params.pop("counter", None)
    if counter:
        copy = [
            {"method": "regexparent", "pattern": prefix}
            for prefix in instrument.counterdict["counters"]
        ]
        target = params.pop("target", None)
        if not target:
            target = f"mean({{{counter}}})"
        expression = f"{{}}*{target}/{{{counter}}}"
        params.update(commonparams)
        ensure_parameter(params, "outputparent", default=coutputparent)
        ensure_parameter(params, "name", "postnormalize")
        task = create_task(
            dependencies=task,
            method="expression",
            expression=expression,
            copy=copy,
            **params,
        )
        tasks.append(task)

    # Remove NaN's
    params = task_parameters(parameters, "replacenan")
    replacenan = params.pop("replacenan", False)
    if replacenan:
        params.update(commonparams)
        ensure_parameter(params, "outputparent", default=coutputparent)
        ensure_parameter(params, "name", "replace")
        ensure_parameter(params, "org", np.nan)
        ensure_parameter(params, "new", 0)
        tmp = create_task(dependencies=task, method="replace", **params)
        tasks.append(tmp)

    # Crop
    params = task_parameters(parameters, "crop")
    cropafter = params.pop("crop", False)
    if cropafter:
        params.update(commonparams)
        ensure_parameter(params, "outputparent", default=coutputparent)
        ensure_parameter(params, "name", "crop")
        ensure_parameter(params, "cropfull", True)
        params["nanfull"] = params.pop("cropfull")
        ensure_parameter(params, "nanval", np.nan)
        ensure_parameter(params, "reference", reference)
        tmp = create_task(dependencies=task, method="crop", **params)
        tasks.append(tmp)

    return tasks

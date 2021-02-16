# -*- coding: utf-8 -*-

import os
from . import batch
from .run import run_sequential
from .id21_ffxas import tasks as ffxas_tasks
from ..utils import instance
from ..instruments.configuration import getinstrument


def staticscan(samplename, datasetname, radix, **parameters):
    jobname = batch.jobname("static", (samplename, datasetname, radix), parameters)

    instrument = getinstrument(**parameters)
    mradix, subdir = instrument.fflocation(samplename, datasetname, type="static")
    if not instance.isarray(radix):
        radix = [radix]
    sourcepaths = [
        os.path.join(parameters["proposaldir"], subdir, rdx) for rdx in radix
    ]

    if len(radix) > 1:
        nxentry = "{}_{}".format(mradix, radix[0], radix[-1])
    else:
        nxentry = "{}_{}".format(mradix, radix[0])
    nxentry = (
        os.path.join(parameters.get("resultsdir", ""), samplename + ".h5")
        + "::/"
        + nxentry
    )
    processdata(jobname, sourcepaths, radix, nxentry, **parameters)


def processdata(jobname, *args, **kwargs):
    if "jobs" in kwargs:
        kwargs["jobs"].append((jobname, processdata_exec, args, kwargs))
    else:
        processdata_exec(*args, **kwargs)


def processdata_exec(sourcepaths, radix, nxentry, **kwargs):
    parameters = dict(kwargs)
    parameters["sourcepaths"] = sourcepaths
    parameters["radix"] = radix
    parameters["nxentry"] = nxentry
    tasks = ffxas_tasks(**parameters)
    if run_sequential(tasks, name="fullfield"):
        pass
    else:
        unfinished = [task for task in tasks if not task.done]
        raise RuntimeError(
            "The following tasks are not finished: {}".format(unfinished)
        )

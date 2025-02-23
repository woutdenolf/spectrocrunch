import os
from ..utils import instance
from ..process.utils import create_task
from ..io import nxfs


def ffparameters(**parameters):
    sourcepaths = parameters["sourcepaths"]
    radix = parameters["radix"]

    rebin = parameters.get("rebin", (1, 1))
    roiraw = parameters.get("roiraw", None)
    stackdim = parameters.get("stackdim", 0)
    normalize = parameters.get("normalize", True)
    normalizeonload = parameters.get("normalizeonload", True)
    flatbeforeafter = parameters.get("flatbeforeafter", True)
    nxentry = parameters.get("nxentry", None)

    if not instance.isarray(sourcepaths):
        sourcepaths = [sourcepaths]
    if not instance.isarray(radix):
        radix = [radix]

    outputparent = nxfs.factory(str(nxentry))
    outputparent = outputparent.parent[outputparent.name + ".1"]

    config = {
        # EDF header
        "frametimelabel": "exposure_time",
        "frametimedefault": 0.0,
        "roilabel": "roi",
        "nbframeslabel": "nb_frames",
        "stacklabel": "energy",
        # Defaults
        "dtype": "np.float32",  # must be floating point!
        "darkcurrentzero": 90.0,
        "darkcurrentgain": 1.0,
        "stackdim": stackdim,
        # Data
        "darklist": list(
            map(
                lambda xy: os.path.join(xy[0], xy[1] + "*_dark_*.edf"),
                zip(sourcepaths, radix),
            )
        ),
        "datalist": list(
            map(
                lambda xy: os.path.join(xy[0], xy[1] + "*_data_*.edf"),
                zip(sourcepaths, radix),
            )
        ),
        "flatlist": list(
            map(
                lambda xy: os.path.join(xy[0], xy[1] + "*_ref_*.edf"),
                zip(sourcepaths, radix),
            )
        ),
        # split up flat fields in before and after (same number of images)
        "flatbeforeafter": flatbeforeafter,
        "normalize": normalize and normalizeonload,
        "keepflat": normalize and not normalizeonload,
        "roi": roiraw,
        "rebin": rebin,
        # Output
        "outputparent": outputparent,
    }

    return config


def tasks(**parameters):
    tasks = []

    # Common parameters
    parameters["stackdim"] = parameters.get("stackdim", 0)
    parameters["default"] = parameters.get("default", "sample")
    commonparams = {k: parameters[k] for k in ["default", "stackdim"]}

    # Image stacks (fullfield)
    ffparams = ffparameters(**parameters)
    ffparams.update(commonparams)
    task = create_task(method="fullfield", name="process:fullfield", **ffparams)
    tasks.append(task)

    # Normalization
    normalize = parameters.get("normalize", False)
    if normalize and not ffparams["normalize"]:
        skip = [
            {"method": "regex", "pattern": "flat1"},
            {"method": "regex", "pattern": "flat2"},
        ]
        if ffparams["flatbeforeafter"]:
            expression = "-ln(2*{}/({flat1}+{flat2}))"
        else:
            expression = "-ln({}/{flat1})"
        task = create_task(
            dependencies=task,
            method="expression",
            name="process:normalize",
            expression=expression,
            skip=skip,
            **commonparams,
        )
        tasks.append(task)

    # Alignment
    alignmethod = parameters.get("alignmethod", None)
    alignreference = "sample"
    if alignmethod and alignreference is not None:
        refimageindex = parameters.get("refimageindex", -1)
        roi = parameters.get("roialign", None)
        plot = parameters.get("plot", False)
        task = create_task(
            dependencies=task,
            method="align",
            name="process:align",
            alignmethod=alignmethod,
            reference=alignreference,
            refimageindex=refimageindex,
            crop=False,
            roi=roi,
            plot=plot,
            **commonparams,
        )
        tasks.append(task)

    # Crop
    roiresult = parameters.get("roiresult", None)
    if roiresult:
        tmp = create_task(
            dependencies=task,
            method="crop",
            name="process:roi",
            roi=roiresult,
            reference=alignreference,
            **commonparams,
        )
        tasks.append(tmp)

    return tasks

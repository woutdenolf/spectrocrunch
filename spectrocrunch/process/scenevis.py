# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import logging
import os

from . import nxprocess
from ..visualization import scene
from ..visualization import scene_view
from ..visualization import defaultsettings
from ..utils import instance
from ..utils import units
from ..io import xiaedf

logger = logging.getLogger(__name__)


def scene_modifyaxes(sc, parameters):
    if parameters.get("originzero", False):
        sc.originzero()
    if parameters.get("transpose", False):
        sc.transpose(True)
    if parameters.get("flipx", False):
        sc.flipx(increasing=parameters.get("incx", None))
    elif parameters.get("incx", False):
        sc.increasingx(True)
    if parameters.get("flipy", False):
        sc.flipy(increasing=parameters.get("incy", None))
    elif parameters.get("incy", False):
        sc.increasingy(True)
    if parameters.get("aspectratio", None) is not None:
        sc.aspectratio = parameters["aspectratio"]


def scene_globalscaling(sc, parameters):
    rlo = parameters.get("rlo", None)
    rhi = parameters.get("rhi", None)
    if rlo is not None or rhi is not None:
        vmin, vmax = sc.vminmax
        dv = np.array(vmax) - np.array(vmin)

        if rlo is None:
            vmin = None
        else:
            rlo, func = instance.asarrayf(rlo)
            vmin = vmin + dv * func(rlo)
        if rhi is None:
            vmax = None
        else:
            rhi, func = instance.asarrayf(rhi)
            vmax = vmin + dv * func(rhi)
        sc.scale(vmin=vmin, vmax=vmax)

    alo = parameters.get("alo", None)
    ahi = parameters.get("ahi", None)
    if alo is not None or ahi is not None:
        sc.scale(vmin=alo, vmax=ahi)


def createwriter(parameters, filename):
    filename, ext = os.path.splitext(filename)
    filename = filename + ".xlsx"
    if filename not in parameters["writers"]:
        writer = pd.ExcelWriter(filename)
        parameters["writers"][filename] = writer
    return parameters["writers"][filename]


def savewriters(parameters):
    for filename, writer in parameters["writers"].items():
        writer.save()
        logger.info("Saved {}".format(filename))


def savefigures(filename, parameters):
    figs = [plt.figure(i) for i in plt.get_fignums()]
    if len(figs) > 1:
        filename, ext = os.path.splitext(filename)
        filename = "{}_{{:02d}}{}".format(filename, ext)
    for i, fig in enumerate(figs):
        name = filename.format(i + 1)
        if not os.path.exists(os.path.dirname(name)):
            os.makedirs(os.path.dirname(name))
        saveparams = parameters.get("saveparams", {})
        fig.savefig(name, **saveparams)
        logger.info("Saved {}".format(name))


def closefigures():
    for i in plt.get_fignums():
        plt.close(i)


def specname(info):
    filename = "spec{:02d}".format(info["spec"])
    radix = "_".join((info["sample"], filename))
    filename = os.path.join(info["sample"], radix, "{}.dat".format(radix))
    return os.path.join(info["root"], filename)


def extractname(item):
    if instance.isarray(item):
        return item[0]
    else:
        return item


def ctrfilenames(info):
    radix = "_".join((info["sample"], info["dataset"]))
    root = info["root"]
    sample = info["sample"]
    fformat = xiaedf.xiaformat_ctr(radix, info["num"])
    return [
        os.path.join(root, sample, radix, "zap", fformat.format(extractname(item)))
        for item in info["items"]
    ]


def specfilename(info):
    filename = "spec{:02d}".format(info["spec"])
    radix = "_".join((info["sample"], filename))
    filename = os.path.join(info["sample"], radix, "{}.dat".format(radix))
    return os.path.join(info["root"], filename)


def createscene(parameters, dependencies=None, output=None):
    parameters = parameters.copy()
    objects = parameters.pop("objects")
    if dependencies is None:
        dependencies = []

    # Create scene and connect it to axes
    figsize = parameters.get("figsize", None)
    if figsize is None:
        figsize = defaultsettings.default("figure.figsize")
    defaultsettings.adapttofigsize(figsize, **parameters)
    f, ax = plt.subplots(figsize=figsize)
    if parameters.get("noaxes", False):
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        # ax.set_axis_off()
    sc = scene.Scene(
        unit0=parameters.get("unitfast", "mm"),
        unit1=parameters.get("unitslow", "mm"),
        title=parameters.get("title", None),
    )
    sc.setaxes(ax)

    # Add objects to scene
    dependency_ctr = -1
    for info in objects:
        plotparams = info.get("plotparams", {}).copy()
        plotparams["scene"] = sc
        dataparams = info.get("dataparams", {}).copy()
        dataparams["instrument"] = parameters.get("instrument", None)
        pmin = plotparams.pop("lo", [])
        pmax = plotparams.pop("hi", [])
        if "spec" in info:
            filename = specname(info)
            item = scene_view.XanesSpec(
                filename, info.get("items", None), plotparams=plotparams, **dataparams
            )
        elif "sample" in info:
            filenames = ctrfilenames(info)
            item = scene_view.ZapRoiMap(
                filenames, info.get("items", None), plotparams=plotparams, **dataparams
            )
        else:
            dependency_ctr += 1
            i = info.get("dependency", dependency_ctr)
            uri = str(dependencies[i])
            item = scene_view.Nexus(
                uri, info.get("items", None), plotparams=plotparams, **dataparams
            )
        if pmin and pmax:
            item.selfscale(pmin=pmin, pmax=pmax)
        item.useaxesnames()

    # Modify axes
    scene_modifyaxes(sc, parameters)

    # Global intensity scaling
    scene_globalscaling(sc, parameters)

    # Save interpolated data
    for item in sc:
        try:
            item.interpolatesave()
        except AttributeError:
            pass

    # Plot/save figure
    sc.updateview()
    if parameters.get("tight_layout", False):
        plt.tight_layout()
    if output is not None:
        savefigures(output, parameters)

    if parameters.get("plot", True):
        # TODO: this hangs when another window is opened in the mean time:
        # sc.ax.get_figure().canvas.mpl_connect('resize_event', lambda event: sc.updateview())
        plt.show()

    closefigures()


class Task(nxprocess.Task):
    """Create scene image
    """

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.required_parameters |= {
            "objects",
            "instrument",
            "figsize",
            "noaxes",
            "unitfast",
            "unitslow",
            "tight_layout",
            "plot",
            "originzero",
            "transpose",
            "flipx",
            "flipy",
            "incx",
            "incy",
            "aspectratio",
            "rlo",
            "rhi",
            "alo",
            "ahi",
            "title",
            "saveparams",
        }
        parameters = self.parameters
        parameters["figsize"] = parameters.get("figsize", None)
        parameters["noaxes"] = parameters.get("noaxes", False)
        parameters["tight_layout"] = parameters.get("tight_layout", False)
        parameters["plot"] = parameters.get("plot", False)
        parameters["unitfast"] = parameters.get("unitfast", "mm")
        parameters["unitslow"] = parameters.get("unitslow", "mm")
        parameters["aspectratio"] = parameters.get("aspectratio", None)
        parameters["originzero"] = parameters.get("originzero", False)
        parameters["transpose"] = parameters.get("transpose", False)
        parameters["flipx"] = parameters.get("flipx", False)
        parameters["flipy"] = parameters.get("flipy", False)
        parameters["incx"] = parameters.get("incx", False)
        parameters["incy"] = parameters.get("incy", False)
        parameters["rlo"] = parameters.get("rlo", None)
        parameters["rhi"] = parameters.get("rhi", None)
        parameters["alo"] = parameters.get("alo", None)
        parameters["ahi"] = parameters.get("ahi", None)
        parameters["title"] = parameters.get("title", None)
        parameters["saveparams"] = parameters.get("saveparams", {})

    def _execute(self):
        createscene(
            self.parameters,
            output=self.temp_localpath.path,
            dependencies=self.previous_outputs,
        )

    @property
    def temp_localpath(self):
        path = super(Task, self).temp_localpath
        return path.parent[path.name + ".png"]

    @property
    def output_localpath(self):
        path = super(Task, self).output_localpath
        return path.parent[path.name + ".png"]

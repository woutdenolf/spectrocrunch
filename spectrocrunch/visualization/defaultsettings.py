import matplotlib.pyplot as plt
import matplotlib
import logging

from ..utils import instance

logger = logging.getLogger(__name__)


def update(parameters):
    matplotlib.rcParams.update(parameters)


def default(name):
    return matplotlib.rcParams[name]


def figsize(publish="screen", aspect=0.75, nsidebyside=1, space=0.0, widescreen=True):
    """
    Args:
        publish(Optional(str)): word, powerpoint, screen
        aspect(Optional(num)): height/width
        nsidebyside(Optional(num)): number of images that should be put side-by-side
        space(Optional(num)): fraction of the total width that is white space
        widescreen(Optional(bool)): only for powerpoint
    Returns:
        tuple: width, height (inch)
    """
    if publish == "word":
        width = 6.24
    elif publish == "powerpoint":
        if widescreen:
            width = 33.867
        else:
            width = 13.33
    else:
        if instance.isnumber(publish):
            width = publish
        else:
            width = 6.4
    width *= (1 - space) / nsidebyside
    return width, width * aspect


def adapttofigsize(size, fontsize=None, fontsizepts=10, **kwargs):
    """
    Args:
        size(2-tuple): inch
        fontsize(Optional(num)): inch
    """
    if fontsize:  # in inch
        # 1 inch = 72 points
        fontsize = fontsize * 72
    else:
        fontsize = min(size) * fontsizepts / 4.8
    linewidth = fontsize * 0.08

    wtick = linewidth
    cwtick = 0.75
    stick = linewidth * 4.375
    cstick = 1 / 1.75
    ptick = linewidth * 4.375
    cptick = 3.4 / 3.5
    parameters = {
        "font.size": fontsize,
        "axes.linewidth": linewidth,
        "lines.linewidth": linewidth * 1.875,
        "patch.linewidth": linewidth * 1.25,
        "axes.titlepad": fontsize * 0.6,
        "axes.labelpad": fontsize * 0.4,
        "xtick.major.width": wtick,
        "xtick.minor.width": wtick * cwtick,
        "ytick.major.width": wtick,
        "ytick.minor.width": wtick * cwtick,
        "xtick.major.size": stick,
        "xtick.minor.size": stick * cstick,
        "ytick.major.size": stick,
        "ytick.minor.size": stick * cstick,
        "xtick.major.pad": ptick,
        "xtick.minor.pad": ptick * cptick,
        "ytick.major.pad": ptick,
        "ytick.minor.pad": ptick * cptick,
        "xtick.direction": kwargs.get("tick.direction", "out"),
        "ytick.direction": kwargs.get("tick.direction", "out"),
    }
    update(parameters)
    parameters = {k: v for k, v in kwargs.items() if k in matplotlib.rcParams}
    update(parameters)


def screensize():
    figs1 = plt.get_fignums()
    mgr = plt.get_current_fig_manager()
    h = mgr.window.winfo_screenheight()
    w = mgr.window.winfo_screenwidth()
    figs2 = plt.get_fignums()
    if len(figs2) > len(figs1):
        plt.close(figs2[-1])
    return w, h


def dpi(publish="photo&text", best=True):
    """
    Args:
        publish(str,tuple): publication medium (WxH when a tuple)
        best(bool):
    """
    if publish == "powerpoint":
        publish = screensize()
    if publish == "color":
        ret = 300
    elif publish == "b&w":
        ret = 600 if best else 300
    elif publish == "photo&text":
        ret = 900 if best else 600
    elif publish == "lineart":
        ret = 1200 if best else 900
    elif instance.issequence(publish):
        width, height = publish
        publish = "width x height = {} x {} inch".format(width, height)
        # Powerpoint: normal and widescreen
        # 50	 500 × 375   667 × 375   50 dpi
        # 96  960 × 720   1280 × 720  96 dpi
        # 100 1000 × 750  1333 × 750  100 dpi
        # 150 1500 × 1125 2000 × 1125 150 dpi
        # 200 2000 × 1500 2667 × 1500 200 dpi
        # 250 2500 × 1875 3333 × 1875 250 dpi
        # 300 3000 × 2250 4000 × 2250 300 dpi
        if best:
            warr = [667, 1280, 133, 200, 2667, 333, 4000]
        else:
            warr = [500, 960, 1000, 1500, 2000, 2500, 3000]
        harr = [375, 720, 750, 1125, 1500, 1875, 2250]
        dpiarr = [50, 96, 100, 150, 200, 250, 300]
        for w, h, r in zip(warr, harr, dpiarr):
            if w > width and h > height:
                ret = r
                break
        else:
            ret = 300
    else:  # screen
        ret = 96
    logger.info("Image publication: {}, {} DPI".format(publish, ret))
    return ret

from collections import OrderedDict
import logging
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as pltcolors
import matplotlib.patches as pltpatches
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ..utils import listtools
from ..utils import instance
from ..utils import units
from ..utils.hashable import CompHashable
from ..math.utils import floatformat
from ..math.utils import round_sig
from ..patch.pint import ureg

logger = logging.getLogger(__name__)


def ColorNorm(name, *args):
    try:
        if name is None:
            name = ""
        name = name.lower()
    except AttributeError:
        return name  # already a normalizer

    if name == "log":

        def norm(**kwargs):
            return pltcolors.LogNorm(**kwargs)

    elif name == "power":

        def norm(**kwargs):
            return pltcolors.PowerNorm(*args, **kwargs)

    elif name == "symlog":

        def norm(**kwargs):
            return pltcolors.SymLogNorm(*args, **kwargs)

    else:

        def norm(**kwargs):
            return pltcolors.Normalize(*args, **kwargs)

    return norm


def ColorNormLinear(name):
    return name != "log" and name != "power" and name != "symlog"


class Geometry2D(object):
    def dataaxis(self, displayaxis):
        """
        Args:
            displayaxis(num): display axis (0 or 1)

        Returns:
            dataaxis(num): data axis (0 or 1)
        """
        if instance.isstring(displayaxis):
            displayaxis = int(displayaxis == "x")

        if self.transposed:
            return 1 - displayaxis
        else:
            return displayaxis

    @property
    def dataaxisx(self):
        """Data axis corresponding to the X-axis"""
        return self.dataaxis("x")

    @property
    def dataaxisy(self):
        """Data axis corresponding to the Y-axis"""
        return self.dataaxis("y")


class Scene(CompHashable, Geometry2D):
    """Each scene can have a number of registered items with associalted plot settings."""

    def __init__(
        self, ax=None, unit0=ureg.dimensionless, unit1=ureg.dimensionless, title=""
    ):
        self.ax = ax
        self._items = OrderedDict()
        self.transposed = False
        self.cmap = plt.get_cmap("jet")
        self.title = title
        self._aspectratio_display = None

        # Order: data dimensions
        self.axlabels = ["Dim0", "Dim1"]
        self._increasing = [True, True]
        self.units = [ureg.Unit(unit0), ureg.Unit(unit1)]
        self.dataoffset = [ureg.Quantity(0, unit0), ureg.Quantity(0, unit1)]
        self.datascale = [1, 1]

    @classmethod
    def _str_unit(cls, unit):
        if unit == ureg.dimensionless:
            return ""
        unit = "{:~}".format(unit)
        if unit == "um":
            unit = r"$\mu$m"
        return unit

    @property
    def strunits(self):
        return [self._str_unit(unit) for unit in self.units]

    def q0(self, v):
        return self.dimquantity(v, 0)

    def q1(self, v):
        return self.dimquantity(v, 1)

    def dimquantity(self, v, dataaxis):
        if isinstance(v, ureg.Quantity):
            return v
        else:
            return ureg.Quantity(v, self.units[dataaxis])

    def xmagnitude(self, v):
        return self._magnitude(v, self.dataaxisx)

    def ymagnitude(self, v):
        return self._magnitude(v, self.dataaxisy)

    def _magnitude(self, v, dataaxis):
        unit = self.units[dataaxis]
        return units.magnitude(v, unit)

    @property
    def xunit(self):
        return self.units[self.dataaxisx]

    @xunit.setter
    def xunit(self, value):
        self.units[self.dataaxisx] = ureg.Unit(value)

    @property
    def yunit(self):
        return self.units[self.dataaxisy]

    @yunit.setter
    def yunit(self, value):
        self.units[self.dataaxisy] = ureg.Unit(value)

    @property
    def strxunit(self):
        return self.strunits[self.dataaxisx]

    @property
    def stryunit(self):
        return self.strunits[self.dataaxisx]

    @property
    def xlabel(self):
        return self._label(self.dataaxisx)

    @property
    def ylabel(self):
        return self._label(self.dataaxisy)

    def _label(self, dataaxis):
        unit = self.strunits[dataaxis]
        label = self.axlabels[dataaxis]
        if unit and label:
            label = "{} ({})".format(label, unit)
        return label

    def setaxes(self, ax):
        self.ax = ax

        # TODO: wsize is already set when callback is called, but asize is not (and we need asize)
        def func(event):
            return self.updateview()

        # def func(event):
        #    print self.wsize_pixels
        #    print self.asize_pixels
        #    self.updateview()

        self.ax.figure.canvas.mpl_connect("resize_event", func)

        # TODO: for scalebar
        # self.ax.callbacks.connect('xlim_changed', func)
        # self.ax.callbacks.connect('ylim_changed', func)

    def showaxes(self):
        self.ax.set_axis_on()
        self.ax.get_dataaxisx().set_visible(True)
        self.ax.get_dataaxisy().set_visible(True)

    def hideaxes(self):
        self.ax.set_axis_off()
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)

    def transpose(self, tr):
        self.transposed = tr

    def register(self, item):
        item.addscene(self)
        self._items[item] = item.defaultsettings()

    def getitemsettings(self, item):
        if item in self._items:
            return self._items[item]
        else:
            return {}

    def __iter__(self):
        for item in self._items:
            item.selectscene(self)
        return iter(self._items)

    @property
    def items(self):
        return list(self)

    def datarange(self, dataaxis):
        ran = [r for item in self for r in item.datarange(dataaxis)]
        return min(ran), max(ran)

    def setdatarange(self, dataaxis, ran):
        self.dataoffset[dataaxis] = self.dimquantity(0, dataaxis)
        self.datascale[dataaxis] = 1
        ranorg = sorted(self.datarange(dataaxis))
        scale = (ran[1] - ran[0]) / (ranorg[1] - ranorg[0])
        off = scale * ranorg[0] - ran[0]
        self.dataoffset[dataaxis] = off
        self.datascale[dataaxis] = scale

    def resetdatarange(self):
        self.dataoffset = [self.q0(0), self.q1(1)]
        self.datascale = [1, 1]

    def increasingx(self, increasing):
        self._increasing[self.dataaxisx] = increasing

    def increasingy(self, increasing):
        self._increasing[self.dataaxisy] = increasing

    def flipx(self, increasing=None):
        self._flipaxis(self.dataaxisx, increasing=increasing)

    def flipy(self, increasing=None):
        self._flipaxis(self.dataaxisy, increasing=increasing)

    def _flipaxis(self, dataaxis, increasing=None):
        if increasing is None:
            self._increasing[dataaxis] = not self._increasing[dataaxis]
        else:
            if increasing == self._increasing[dataaxis]:
                self.setdatarange(dataaxis, self.datarange(dataaxis)[::-1])
            else:
                self._increasing[dataaxis] = not self._increasing[dataaxis]

    def origin(self, dataaxis, off):
        self.dataoffset[dataaxis] = off - self.dataoffset[dataaxis]

    def originreset(self):
        self.origin(0, self.q0(0))
        self.origin(1, self.q1(0))

    def originzero(self):
        self.origin(0, self.datarange(0)[0])
        self.origin(1, self.datarange(1)[0])

    def datatransform(self, arr, dataaxis):
        arr = units.asqarray(arr)
        return arr * self.datascale[dataaxis] - self.dataoffset[dataaxis]

    def displaylim(self, lim, dataaxis):
        if self._increasing[dataaxis]:
            return min(lim), max(lim)
        else:
            return max(lim), min(lim)

    @property
    def displaylimx(self):
        datalimx = [x for item in self for x in item.datalimx]
        return self.displaylim(datalimx, self.dataaxisx)

    @property
    def displaylimy(self):
        ylim = [y for item in self for y in item.datalimy]
        return self.displaylim(ylim, self.dataaxisy)

    @property
    def aspectcorrect(self):
        if self._aspectratio_display is None:
            return None  # same as 1: aspect ratio of display and data are the same
        return self._aspectratio_display / self.aspectratio_data

    @property
    def aspectratio(self):
        if self._aspectratio_display is None:
            return self.aspectratio_data
        else:
            return self._aspectratio_display

    @aspectratio.setter
    def aspectratio(self, value):
        self._aspectratio_display = value

    @property
    def aspectratio_data(self):
        displaylimy = self.displaylimy
        displaylimx = self.displaylimx
        dx = abs(displaylimx[1] - displaylimx[0])
        dy = abs(displaylimy[1] - displaylimy[0])
        return (dy / dx).to("dimensionless").magnitude

    def crop(self):
        lim = self.xmagnitude(self.displaylimx)
        self.ax.set_xlim(lim[0], lim[1])
        # self.ax.set_xbound(lim[0],lim[1])
        lim = self.ymagnitude(self.displaylimy)
        self.ax.set_ylim(lim[0], lim[1])
        # self.ax.set_ybound(lim[0],lim[1])

    def updateview(self, **kwargs):
        for item in self:
            item.updateview(**kwargs)
        self.crop()

        if self.title is not None:
            self.ax.set_title(self.title)
        if self.xlabel is not None:
            self.ax.set_xlabel(self.xlabel)
        if self.ylabel is not None:
            self.ax.set_ylabel(self.ylabel)

    @property
    def vminmax(self):
        vmin = None
        vmax = None
        tmp = [item.vminmax for item in self if hasattr(item, "vminmax")]
        if not tmp:
            return None, None

        vmin, vmax = zip(*tmp)
        vmin = [listtools.aslist(v) for v in vmin]
        vmax = [listtools.aslist(v) for v in vmax]
        n = max(len(v) for v in vmin)

        def reshape(v):
            if len(v) != n:
                v.extend([np.nan] * (n - len(v)))
            return v

        vmin = np.array([reshape(v) for v in vmin])
        vmax = np.array([reshape(v) for v in vmax])
        vmin = np.nanmin(vmin, axis=0)
        vmax = np.nanmax(vmax, axis=0)
        return vmin, vmax

    def scale(self, vmin=None, vmax=None):
        self.set_setting("vmin", vmin)
        self.set_setting("vmax", vmax)

    def set_setting(self, attr, value):
        for item, settings in self._items.items():
            if attr in settings:
                settings[attr] = value

    def set_settings(self, settings):
        for item, _settings in self._items.items():
            _settings.update(settings)

    @property
    def wsize_pixels(self):
        bbox = self.ax.figure.get_window_extent().transformed(
            self.ax.figure.dpi_scale_trans.inverted()
        )
        width, height = (
            bbox.width * self.ax.figure.dpi,
            bbox.height * self.ax.figure.dpi,
        )
        return width, height

    @property
    def asize_pixels(self):
        bbox = self.ax.get_window_extent().transformed(
            self.ax.figure.dpi_scale_trans.inverted()
        )
        width, height = (
            bbox.width * self.ax.figure.dpi,
            bbox.height * self.ax.figure.dpi,
        )
        return width, height

    @property
    def asize_data(self):
        x1, x2 = self.xmagnitude(self.displaylimx)
        y1, y2 = self.ymagnitude(self.displaylimy)
        return x2 - x1, y2 - y1

    def fontsize_data(self, pt):
        xp, yp = self.asize_pixels
        xd, yd = self.asize_data
        return pt * xd / float(xp), pt * yd / float(yp)


class Item(CompHashable, Geometry2D):
    """Each item can be registered with multiple scenes. Data is owned by the Item instance,
    plot settings are owned by the scenes.
    """

    def __init__(self, scene=None, name=None, **kwargs):
        self._scenes = OrderedDict()
        self._sceneindex = -1
        self.name = name

        if scene is not None:
            self.register(scene)

        Item._updatedata(self, **kwargs)

    def _updatedata(self, axis0name="Dim0", axis1name="Dim1", **settings):
        self.axis0name = axis0name
        self.axis1name = axis1name
        try:
            if settings:
                self.set_settings(settings)
        except RuntimeError:
            logger.warning("Item settings are not applied (provide a scene)")

    def updatedata(self, **kwargs):
        Item._updatedata(self, **kwargs)

    def register(self, scene):
        scene.register(self)

    def useaxesnames(self):
        self.scene.axlabels = [self.axis0name, self.axis1name]

    def defaultsettings(self):
        return {}

    def set_settings(self, settings):
        settings2 = self.scene.getitemsettings(self)
        for k, v in settings.items():
            settings2[k] = v

    def get_settings(self, keys):
        settings = self.scene.getitemsettings(self)
        return {k: settings[k] for k in keys}

    def set_setting(self, attr, value):
        settings = self.scene.getitemsettings(self)
        settings[attr] = value

    def get_setting(self, key):
        settings = self.scene.getitemsettings(self)
        return settings[key]

    def selectscene(self, s):
        try:
            self._sceneindex = list(self._scenes.keys()).index(s)
        except Exception:
            raise RuntimeError("This object is not registered with scene {}".format(s))

    @property
    def scene(self):
        if len(self._scenes) == 0:
            raise RuntimeError("This object is not registered with any scene")
        return list(self._scenes.keys())[self._sceneindex]

    def addscene(self, s):
        if s not in self._scenes:
            self._scenes[s] = OrderedDict()

    @property
    def sceneitems(self):
        """My items in the active scene"""
        return self._scenes.get(self.scene, OrderedDict())

    def removefromscene(self):
        """Remove myself from the active scene"""
        items = self.sceneitems
        for item in items:
            if item is not None:
                items[item].remove()
        self._scenes[self.scene] = OrderedDict()

    def refreshscene(self, newitems):
        """Update the active scene with new items"""
        olditems = self.sceneitems

        for name in olditems:
            if name in newitems:
                if newitems[name] != olditems[name]:  # TODO: use "is not"?
                    olditems[name].remove()
            else:
                olditems[name].remove()
        self._scenes[self.scene] = newitems

    def datarange(self, dataaxis):
        raise NotImplementedError("Item is an abstract class")

    @property
    def datalimx(self):
        return self.datarange(self.dataaxisx)

    @property
    def datalimy(self):
        return self.datarange(self.dataaxisy)

    @property
    def transposed(self):
        return self.scene.transposed


class Image(Item):
    def __init__(self, data, lim0=None, lim1=None, labels=None, **kwargs):
        super(Image, self).__init__(**kwargs)
        Image._updatedata(self, data, lim0=lim0, lim1=lim1, labels=labels)

    def _updatedata(self, data, lim0=None, lim1=None, labels=None):
        self.data = data
        self.labels = labels

        if lim0 is None:
            lim0 = [0, data.shape[0] - 1]
        else:
            lim0 = [lim0[0], lim0[-1]]

        if lim1 is None:
            lim1 = [0, data.shape[1] - 1]
        else:
            lim1 = [lim1[0], lim1[-1]]

        self.lim = [lim0, lim1]

    def updatedata(self, data, lim0=None, lim1=None, labels=None, **settings):
        Image._updatedata(self, data, lim0=lim0, lim1=lim1, labels=labels)
        super(Image, self).updatedata(**settings)

    @property
    def data(self):
        """Data for the selected channels"""
        channels = self.channels
        shape = (self._data.shape[0], self._data.shape[1], len(channels))
        data = np.zeros(shape, self.dtype)
        for ch, i in enumerate(channels):
            if i is not None:
                data[..., ch] = self._get_data(i)
        if len(channels) == 1:
            data = data[..., 0]
        u = self.get_setting("colorbar_units")
        if u:
            data = data / float(units.Quantity(1, u).to("dimensionless").magnitude)
        return data

    def _get_data(self, index):
        if self._data.ndim == 3:
            return self._data[..., index]
        else:
            if index != 0:
                raise IndexError(
                    "Raw data has only one channel, you tried to select channel {}".format(
                        index
                    )
                )
            return self._data

    @data.setter
    def data(self, value):
        self._data = value

    @property
    def dtype(self):
        return self._data.dtype

    @property
    def labels(self):
        """Labels for the selected channels"""
        labels = [None] * self.nchannels
        for ch, i in enumerate(self.channels):
            if i is not None:
                try:
                    labels[ch] = self._labels[i]
                except Exception:
                    labels[ch] = "<empty>"
        return labels

    @labels.setter
    def labels(self, value):
        self._labels = value

    @property
    def channels(self):
        """Return 1 or 3 channels"""
        settings = self.scene.getitemsettings(self)
        channels = settings["channels"]
        if channels is None:
            nchannels = min(self.datadepth, 3)
            channels = range(nchannels)
        else:
            if not instance.isarray(channels):
                channels = [channels]
            nchannels = min(len(channels), 3)
        if nchannels != 1:
            nchannels = 3

        if len(channels) < nchannels:
            channels.extend([None] * (nchannels - len(channels)))
        else:
            channels = channels[0:nchannels]
        return channels

    @property
    def datadepth(self):
        if self._data.ndim == 3:
            return self._data.shape[-1]
        else:
            return 1

    @property
    def nchannels(self):
        return len(self.channels)

    @property
    def datadisplay(self):
        if self.transposed:
            return np.swapaxes(self.data, 0, 1)
        else:
            return self.data

    @property
    def vminmax(self):
        ptrim = max(min(self.get_setting("ptrim"), 1), 0) * 100
        img = self.data
        if self.data.ndim == 3:
            nchannels = self.data.shape[-1]
            vmin = np.asarray(
                [
                    np.nanpercentile(img[..., i], ptrim, interpolation="nearest")
                    for i in range(nchannels)
                ]
            )
            vmax = np.asarray(
                [
                    np.nanpercentile(img[..., i], 100 - ptrim, interpolation="nearest")
                    for i in range(nchannels)
                ]
            )
        else:
            vmin = np.nanpercentile(img, ptrim, interpolation="nearest")
            vmax = np.nanpercentile(img, 100 - ptrim, interpolation="nearest")
        return vmin, vmax

    def scale(self, vmin=None, vmax=None):
        self.set_setting("vmin", vmin)
        self.set_setting("vmax", vmax)

    @staticmethod
    def _pminmaxparse(p, v):
        # Handle mismatch in number of channels
        p = instance.asscalar(p)
        if instance.isarray(p):
            p = instance.asarray(p)
            if instance.isarray(v):
                if p.size < v.size:
                    p = np.append(p, [0] * (v.size - p.size))
                else:
                    p = p[0 : v.size]
            else:
                p = p[0]
        return p, v

    def selfscale(self, pmin=0, pmax=1, ptrim=0):
        self.set_setting("ptrim", ptrim)
        vmin, vmax = self.vminmax
        d = vmax - vmin
        pmin, vmin = self._pminmaxparse(pmin, vmin)
        pmax, vmax = self._pminmaxparse(pmax, vmax)
        vmax = vmin + d * pmax
        vmin = vmin + d * pmin
        self.scale(vmin=vmin, vmax=vmax)

    @staticmethod
    def defaultsettings():
        return {
            "cmap": None,
            "aspectcorrect": None,
            "border": False,
            "color": None,
            "image": True,
            "linewidth": 2,
            "alpha": 1,
            "vmin": None,
            "vmax": None,
            "ptrim": 0.0,
            "cnorm": None,
            "cnormargs": (),
            "datafunc": None,
            "legend": False,
            "legend_position": "RT",
            "legend_spacing": 2.0,
            "scalebar": False,
            "scalebar_position": "lower right",
            "scalebar_ratio": 1 / 5.0,
            "scalebar_fraction": 0.1,
            "scalebar_pad": 0.1,
            "scalebar_labelpad": 2,
            "scalebar_color": "#ffffff",
            "scalebar_size": 0.0,
            "fontsize": matplotlib.rcParams["font.size"],
            "fontweight": 500,
            "channels": None,
            "colorbar": False,
            "colorbar_units": "",
            "compositions": {},
            "title": None,
        }

    def datarange(self, dataaxis, border=False):
        lim = self.scene.datatransform(self.lim[dataaxis], dataaxis)
        if border:
            d = (lim[1] - lim[0]) / (2.0 * (self._data.shape[dataaxis] - 1))
            return lim[0] - d, lim[1] + d
        else:
            return lim[0], lim[1]

    @property
    def datalimx(self):
        return self.datarange(self.dataaxisx, border=True)

    @property
    def datalimy(self):
        return self.datarange(self.dataaxisy, border=True)

    @property
    def displaylimx(self):
        scene = self.scene
        return scene.displaylim(self.datalimx, scene.dataaxisx)

    @property
    def displaylimy(self):
        scene = self.scene
        return scene.displaylim(self.datalimy, scene.dataaxisy)

    @property
    def labelxy(self):
        scene = self.scene
        settings = scene.getitemsettings(self)

        if "S" in settings["legend_position"]:
            o = scene
        else:
            o = self

        if "L" in settings["legend_position"]:
            x = 0
            horizontalalignment = "right"
            ox = -0.5
            dx = 0
        else:
            x = 1
            horizontalalignment = "left"
            ox = 0.5
            dx = 0

        if "B" in settings["legend_position"]:
            y = 0
            verticalalignment = "bottom"
            oy = 0
            dy = settings["legend_spacing"]
        else:
            y = 1
            verticalalignment = "top"
            oy = 0
            dy = -settings["legend_spacing"]

        x = scene.xmagnitude(o.displaylimx[x])
        y = scene.ymagnitude(o.displaylimy[y])
        xchar, ychar = scene.fontsize_data(settings["fontsize"])

        return (
            x + xchar * ox,
            y + ychar * oy,
            xchar * dx,
            ychar * dy,
            horizontalalignment,
            verticalalignment,
        )

    def bordercoord(self):
        x = sorted(self.datalimx)
        y = sorted(self.datalimy)
        x = [x[0], x[1], x[1], x[0], x[0]]
        y = [y[0], y[0], y[1], y[1], y[0]]
        return x, y

    @classmethod
    def compositioncolor(cls, name, v, vmin, vmax):
        color = [0, 0, 0]
        v2 = [0, 0, 0]
        for i, (c, mi, ma) in enumerate(zip(v, vmin, vmax)):
            v2[i] = c
            if mi is not None and ma is not None:
                color[i] = max(min((c - mi) / (ma - mi), 1), 0)
                v2[i] = color[i] * (ma - mi) + mi
        return name, v2, color

    def updateview(self):
        scene = self.scene

        settings = scene.getitemsettings(self)
        if settings["cmap"] is None:
            cmap = scene.cmap
        else:
            cmap = settings["cmap"]
        if settings["aspectcorrect"] is None:
            aspectcorrect = scene.aspectcorrect
        else:
            aspectcorrect = settings["aspectcorrect"]
        alpha = settings["alpha"]

        # Create intensity-scaled image
        vmin, vmax = scene.vminmax
        if settings["vmin"] is not None:
            vmin = settings["vmin"]
        if settings["vmax"] is not None:
            vmax = settings["vmax"]
        image = self.datadisplay

        # clip = True -> neglect colormap's over/under/masked colors
        nchannels = self.nchannels
        if image.ndim == 3:
            if not instance.isarray(vmin):
                vmin = [vmin] * nchannels
            if not instance.isarray(vmax):
                vmax = [vmax] * nchannels
            normcb = []

            for i in range(nchannels):
                norm = ColorNorm(settings["cnorm"], *settings["cnormargs"])(
                    vmin=vmin[i], vmax=vmax[i], clip=True
                )
                image[..., i] = norm(image[..., i]).data
                normcb.append(norm)
            norm = None
        else:
            if instance.isarray(vmin):
                vmin = vmin[0]
            if instance.isarray(vmax):
                vmax = vmax[0]
            norm = ColorNorm(settings["cnorm"], *settings["cnormargs"])(
                vmin=vmin, vmax=vmax, clip=True
            )
            normcb = []

        # Coordinates and borders
        extent = list(scene.xmagnitude(self.datalimx)) + list(
            scene.ymagnitude(self.datalimy)
        )
        x, y = self.bordercoord()
        x = scene.xmagnitude(x)
        y = scene.ymagnitude(y)

        # Create/update objects
        items = self.sceneitems
        newitems = OrderedDict()

        # Image
        if "image" in items:
            newitems["image"] = items["image"]
            plt.setp(
                newitems["image"],
                interpolation="nearest",
                extent=extent,
                norm=norm,
                alpha=alpha,
                cmap=cmap,
            )
        else:
            newitems["image"] = scene.ax.imshow(
                image,
                interpolation="nearest",
                extent=extent,
                cmap=cmap,
                origin="lower",
                norm=norm,
                alpha=alpha,
                aspect=aspectcorrect,
            )
        newitems["image"].set_visible(settings["image"])

        # Border
        if "border" in items:
            newitems["border"] = items["border"]
            newitems["border"].set_data(x, y)
            kwargs = {k: settings[k] for k in ["linewidth", "color", "alpha"]}
            plt.setp(newitems["border"], **kwargs)
        else:
            kwargs = {k: settings[k] for k in ["linewidth", "color", "alpha"]}
            newitems["border"] = scene.ax.plot(x, y, **kwargs)[0]
        newitems["border"].set_visible(settings["border"])
        settings["color"] = newitems["border"].get_color()

        # Colorbar
        legend, colors, lvisible = self.addcolorbar(
            scene, settings, items, newitems, normcb, vmin, vmax, color=kwargs["color"]
        )

        # Add composition to legend
        self.addcomposition(settings, legend, colors, vmin, vmax)

        # Legend
        self.addlegend(
            scene, settings, items, newitems, legend, colors, visible=lvisible
        )

        # Scalebar
        self.addscalebar(
            scene,
            items,
            newitems,
            position=settings["scalebar_position"],
            visible=settings["scalebar"],
            ratio=settings["scalebar_ratio"],
            xfrac=settings["scalebar_fraction"],
            pad=settings["scalebar_pad"],
            color=settings["scalebar_color"],
            size=settings["scalebar_size"],
            sep=settings["scalebar_labelpad"],
        )

        self.refreshscene(newitems)

    def addcomposition(self, settings, legend, colors, vmin, vmax):
        compositions = []
        if self.nchannels > 1:
            for name, comp in settings["compositions"].items():
                compositions.append(self.compositioncolor(name, comp, vmin, vmax))
        for name, value, color in compositions:
            if name not in legend:
                legend.append(name)
                colors.append(color)

    def addcolorbar(
        self, scene, settings, items, newitems, normcb, vmin, vmax, color=None
    ):
        visible = settings["legend"]

        legend = []

        if settings["colorbar"]:
            # TODO: update existing color bar?
            if normcb:
                # Triangle: cannot handle normcb which is not linear
                # if settings["colorbar"]==2 and ColorNormLinear(settings["cnorm"]):
                #    #This assumes sum is 1 so not useful here
                #    colorbar_rgb.triangle(vmin=vmin,vmax=vmax,names=self.labels,compositions = compositions)
                # else:
                #    newitems.update(colorbar_rgb.bars(vmin=vmin,vmax=vmax,names=legend,norms=normcb,ax=scene.ax))
                #    legend = []
                legend = self.labels
                for i, (name, mi, ma) in enumerate(zip(legend, vmin, vmax)):
                    if name is not None:
                        if ma > 0.01:
                            fmt = floatformat(ma, 3)
                        else:
                            fmt = ":.2E"

                        visible = True
                        u = settings["colorbar_units"]
                        if u:
                            u = "{:~}".format(units.ureg.Unit(u))
                            u = " " + u.replace(" ", "")
                        else:
                            u = ""
                        if settings["legend"]:
                            fmt = "{{}} [{{{}}},{{{}}}]{}".format(fmt, fmt, u)
                            legend[i] = fmt.format(name, mi, ma)
                        else:
                            fmt = "[{{{}}},{{{}}}]{}".format(fmt, fmt, u)
                            legend[i] = fmt.format(mi, ma)
            else:
                divider = make_axes_locatable(scene.ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                newitems["colorbar"] = plt.colorbar(
                    newitems["image"], label=settings["colorbar_units"], cax=cax
                )

        nchannels = self.nchannels
        if nchannels == 1:
            colors = [color]
        else:
            colors = ["r", "g", "b"]
            if nchannels > 3:
                colors.extend([None] * (nchannels - 3))
            colors = colors[0:nchannels]

        if settings["title"]:
            colors = [settings["color"]] + colors
            legend = [settings["title"]] + legend

        return legend, colors, visible

    def addlegend(
        self,
        scene,
        settings,
        items,
        newitems,
        labels,
        colors,
        visible=True,
        astext=False,
        withpatch=False,
    ):
        if not labels:
            return
        astext = len(labels) == 1

        if astext:
            (
                xlabel,
                ylabel,
                dx,
                dy,
                horizontalalignment,
                verticalalignment,
            ) = self.labelxy
            kwargs = {k: settings[k] for k in ["color", "alpha"]}
            kwargs["horizontalalignment"] = horizontalalignment
            kwargs["verticalalignment"] = verticalalignment
            i = -1
        else:
            kwargs = {k: settings[k] for k in ["color"]}
            patches = []
            if not withpatch:
                labels2 = []
                colors2 = []

        for label, color in zip(labels, colors):
            if label is None or color is None:
                continue

            kwargs["color"] = color

            if astext:
                i += 1

                name = "label{}".format(i)
                if name in items:
                    newitems[name] = items[name]
                    newitems[name].set_text(label)
                    newitems[name].set_position((xlabel + i * dx, ylabel + i * dy))
                    plt.setp(newitems[name], **kwargs)
                else:
                    newitems[name] = scene.ax.text(
                        xlabel + i * dx, ylabel + i * dy, label, **kwargs
                    )
                newitems[name].set_visible(visible)
            else:
                if withpatch:
                    p = pltpatches.Patch(label=label, **kwargs)
                else:
                    labels2.append(label)
                    colors2.append(color)
                    p = pltpatches.Rectangle(
                        (0, 0), 1, 1, fill=False, visible=False, **kwargs
                    )
                patches.append(p)

        if not astext and patches:
            bbox_to_anchor = [1, 1]
            # if "R" in settings["legend_position"]:
            #    bbox_to_anchor[0] = 1.05

            if withpatch:
                legend = plt.legend(
                    handles=patches,
                    loc="upper left",
                    frameon=False,
                    borderaxespad=0,
                    bbox_to_anchor=bbox_to_anchor,
                )
            else:
                legend = plt.legend(
                    handles=patches,
                    labels=labels2,
                    loc="upper left",
                    frameon=False,
                    borderaxespad=0,
                    bbox_to_anchor=bbox_to_anchor,
                    handlelength=0,
                    handletextpad=0,
                )
                for text, c in zip(legend.get_texts(), colors2):
                    plt.setp(text, color=c)

            legend.set_visible(visible)
            newitems["legend"] = legend

    def addscalebar(
        self,
        scene,
        items,
        newitems,
        position="lower right",
        visible=True,
        ratio=8.0,
        xfrac=0.1,
        pad=0.1,
        sep=2,
        color="#ffffff",
        size=0,
        unit=None,
    ):
        if not visible:
            return

        name = "scalebar"
        if unit is None:
            unit = scene.xunit

        # Scalebar horizontal size
        sx = scene.xmagnitude(self.displaylimx)
        sx = units.Quantity(abs(sx[1] - sx[0]), scene.xunit)
        if size:
            xsize = units.Quantity(size, scene.xunit)
        else:
            xsize = units.Quantity(round_sig(sx.magnitude * xfrac, 1), scene.xunit)
        xfrac = (xsize / sx).to("dimensionless").magnitude

        # Integer if possible
        xsize = xsize.to(unit)
        mxsize = xsize.magnitude
        if mxsize == int(mxsize):
            mxsize = int(mxsize)
            xsize = units.Quantity(mxsize, xsize.units)

        # Scalebar vertical size in axis units
        ys = scene.ymagnitude(self.displaylimy)
        ys = abs(ys[1] - ys[0])
        yfrac = ratio * xfrac
        mysize = ys * yfrac

        # TODO: handle zoom

        scalebar = AnchoredSizeBar(
            scene.ax.transData,
            mxsize,
            "{:~}".format(xsize),
            position,
            pad=pad,
            color=color,
            sep=sep,
            fill_bar=True,
            frameon=False,
            size_vertical=mysize,
        )

        newitems[name] = scalebar
        newitems[name].set_visible(visible)

        scene.ax.add_artist(scalebar)


class Scatter(Item):
    def __init__(self, coordinate0, coordinate1, labels=None, **kwargs):
        super(Scatter, self).__init__(**kwargs)
        Scatter._updatedata(self, coordinate0, coordinate1, labels=labels)

    def _updatedata(self, coordinate0, coordinate1, labels=None):
        self._coordinate = [
            instance.asarray(coordinate0),
            instance.asarray(coordinate1),
        ]
        # self._coordinate = [units.asqarray(coordinate0),units.asqarray(coordinate1)]
        # self._coordinate = [coordinate0,coordinate1]
        if labels is None:
            self.labels = []
        else:
            self.labels = labels

    def updatedata(self, coordinate0, coordinate1, labels=None, **settings):
        Image._updatedata(self, coordinate0, coordinate1, labels=labels)
        super(Scatter, self).updatedata(**settings)

    def datarange(self, dataaxis):
        ran = self.scene.datatransform(self._coordinate[dataaxis], dataaxis)
        return min(ran), max(ran)

    def coordinate(self, dataaxis):
        return self.scene.datatransform(self._coordinate[dataaxis], dataaxis)

    @staticmethod
    def defaultsettings():
        return {
            "color": None,
            "marker": "+",
            "linestyle": "",
            "linewidth": 2,
            "fill": False,
            "alpha": None,
            "closed": False,
            "labels": True,
            "skiplabels": [],
            "horizontalalignment": "left",
            "verticalalignment": "bottom",
            "fontsize": matplotlib.rcParams["font.size"],
            "labeloffset": 0.1,
            "fontweight": 500,
            "mew": 0.5,
            "ms": 6,
        }

    def updateview(self):
        scene = self.scene
        settings = scene.getitemsettings(self)

        x = instance.asarray(scene.xmagnitude(self.coordinate(self.dataaxisx)))
        y = instance.asarray(scene.ymagnitude(self.coordinate(self.dataaxisy)))

        # Create/update objects
        items = self.sceneitems
        newitems = OrderedDict()

        # Polygon
        if settings["fill"]:
            kwargs = {
                k: settings[k]
                for k in ["linestyle", "linewidth", "color", "alpha", "closed"]
            }
            if not kwargs["linestyle"]:
                kwargs["linestyle"] = None

            if "patch" in items:
                newitems["patch"] = items["patch"]
                newitems["patch"].set_xy = zip(x, y)
                plt.setp(newitems["patch"], **kwargs)
            else:
                newitems["patch"] = matplotlib.patches.Polygon(zip(x, y), **kwargs)
                scene.ax.add_patch(newitems["patch"])
            settings["color"] = newitems["patch"].get_facecolor()

        # Polyline/points
        if settings["linestyle"] or settings["marker"]:
            kwargs = {
                k: settings[k]
                for k in [
                    "marker",
                    "linestyle",
                    "alpha",
                    "color",
                    "linewidth",
                    "mew",
                    "ms",
                ]
            }

            if settings["closed"] and settings["linestyle"]:
                xcl = np.append(x, x[0])
                ycl = np.append(y, y[0])
            else:
                xcl = x
                ycl = y

            if "line" in items:
                newitems["line"] = items["line"]
                newitems["line"].set_data(xcl, ycl)
                plt.setp(newitems["line"], **kwargs)

            else:
                newitems["line"] = scene.ax.plot(xcl, ycl, **kwargs)[0]
            settings["color"] = newitems["line"].get_color()

        # Labels
        if settings["labels"] and bool(self.labels):
            xo, yo = scene.fontsize_data(settings["fontsize"])
            xo *= settings["labeloffset"]
            yo *= settings["labeloffset"]

            kwargs = {
                k: settings[k]
                for k in [
                    "horizontalalignment",
                    "verticalalignment",
                    "alpha",
                    "color",
                    "fontsize",
                    "fontweight",
                ]
            }

            for i, (xi, yi, label) in enumerate(zip(x, y, self.labels)):
                if label in settings["skiplabels"]:
                    label = ""
                name = "label{}".format(i)
                if name in items:
                    newitems[name] = items[name]
                    newitems[name].set_text(label)
                    newitems[name].set_position((xi, yi))
                    newitems[name].xytext = (xi + xo, yi + yo)  # has no effect!!!!
                    plt.setp(newitems[name], **kwargs)
                else:
                    newitems[name] = scene.ax.annotate(
                        label,
                        xy=(xi, yi),
                        xytext=(xi + xo, yi + yo),
                        xycoords="data",
                        textcoords="data",
                        **kwargs,
                    )

        self.refreshscene(newitems)


class Polyline(Scatter):
    @staticmethod
    def defaultsettings():
        settings = Scatter.defaultsettings()
        settings["marker"] = ""
        settings["linestyle"] = "-"
        settings["closed"] = True
        return settings


class Polygon(Scatter):
    @staticmethod
    def defaultsettings():
        settings = Scatter.defaultsettings()
        settings["marker"] = ""
        settings["fill"] = True
        settings["closed"] = True
        return settings


class Text(Scatter):
    @staticmethod
    def defaultsettings():
        settings = Scatter.defaultsettings()
        settings["labels"] = True
        return settings

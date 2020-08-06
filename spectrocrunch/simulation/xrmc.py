# -*- coding: utf-8 -*-

import os
import re
import itertools
import logging
import functools
import numpy as np
from glob import glob
from collections import MutableMapping, OrderedDict
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from future.utils import with_metaclass
from abc import ABCMeta, abstractmethod, abstractproperty
from scipy.spatial.transform import Rotation

from ..io import localfs
from ..io import mca
from ..utils import subprocess
from ..utils import instance
from ..math import quadrics
from ..materials.compoundfromname import compoundfromname


logger = logging.getLogger(__name__)


def installed():
    return subprocess.installed("xrmc")


def execute(*args, **kwargs):
    out, err, returncode = subprocess.execute("xrmc", *args, **kwargs)
    if isinstance(out, bytes):
        out = out.decode()
    if isinstance(err, bytes):
        err = err.decode()
    return out, err, returncode


headerdt = np.dtype(
    [
        ("ninter", np.int32),  # N. of interactions
        ("ncol", np.int32),  # N. of columns
        ("nrow", np.int32),  # N. of rows
        ("scol", np.float64),  # Pixel size Sx
        ("srow", np.float64),  # Pixel size Sy
        ("time", np.float64),  # Exposure time in sec.
        ("type", np.int32),  # Pixel content type
        ("nbin", np.int32),  # N. of energy bins
        ("emin", np.float64),  # Minimum bin energy
        ("emax", np.float64),
    ]
)  # Maximum bin energy


def saveemptyresult(filename, header):
    with open(filename, "w") as f:
        header.tofile(f)
        shape = header["ninter"], header["nbin"], header["ncol"], header["nrow"]
        np.zeros(np.product(shape), dtype=np.float64).tofile(f)


def _loadxrmcresult(filename):
    """
    returns: Ninteractions, Nrows, Ncolumns, Nchannels
    """
    with open(filename, "rb") as f:
        header = np.fromfile(f, dtype=headerdt, count=1)[0]
        shape = header["ninter"], header["nbin"], header["ncol"], header["nrow"]
        data = np.fromfile(f, dtype=np.float64).reshape(shape)
        x0 = np.arange(header["nrow"]) - (header["nrow"] - 1) / 2.0
        x1 = np.arange(header["ncol"]) - (header["ncol"] - 1) / 2.0
        x0 *= header["srow"]
        x1 *= header["scol"]
        energy = np.linspace(header["emin"], header["emax"], header["nbin"])
        zero = header["emin"]
        gain = (header["emax"] - header["emin"]) / (header["nbin"] - 1.0)
        info = {
            "time": header["time"],
            "x0": x0,
            "x1": x1,
            "xenergy": energy,
            "zero": zero,
            "gain": gain,
        }
    return np.transpose(data, (0, 3, 2, 1)), info


def loadxrmcresult(path, basename, ext=".dat"):
    """
    returns: Ninteractions, Nstack, Nrows, Ncolumns, Nchannels
    """
    filenames = glob(os.path.join(path, basename + "*" + ext))
    data = []
    info = {}
    for filename in sorted(filenames):
        datai, infoi = _loadxrmcresult(filename)
        data.append(datai)
        if not info:
            info = infoi
    data = np.stack(data, axis=1)
    return data, info


def xrmcresult_to_mca(xrmcfile, mcafile, mode="w"):
    """
    Save sum spectrum as mca
    """
    data, info = _loadxrmcresult(xrmcfile)
    mcadata = data.sum(axis=tuple(range(data.ndim - 1)))
    mca.save(mcadata, mcafile, mode=mode, zero=info["zero"], gain=info["gain"])


def showxrmcresult(
    data, x0=None, x1=None, xenergy=None, time=None, ylog=False, **kwargs
):
    print("Data shape: {}".format(data.shape))

    if len(x0) > 1:
        x0min, x0max = x0[0], x0[-1]
        dx0 = (x0[1] - x0[0]) * 0.5
        x0min -= dx0
        x0max += dx0
        x1min, x1max = x1[0], x1[-1]
        dx1 = (x1[1] - x1[0]) * 0.5
        x1min -= dx1
        x1max += dx1
        extent = x1min, x1max, x0min, x0max
    else:
        extent = None

    axes = tuple(range(2, data.ndim))
    min = data.min(axes).T
    max = data.max(axes).T
    sum = data.sum(axes).T
    for i, (mi, ma, su) in enumerate(zip(min, max, sum)):
        print("Step {}".format(i))
        print(" Min counts/pixel (for each order): {}".format(mi))
        print(" Max counts/pixel (for each order): {}".format(ma))
        print(" Total counts (for each order): {}".format(su))

    for order, stack in enumerate(data):
        for islice, chunk in enumerate(stack):
            chunk = np.squeeze(chunk)
            if chunk.ndim == 0:
                # Diode signal
                continue
            elif chunk.ndim == 1:
                # XRF spectrum
                spectrum = chunk
                img = None
            elif chunk.ndim == 2:
                # Transmission image
                spectrum = None
                img = chunk
            else:
                # XRF image
                img = chunk.sum(axis=-1)
                spectrum = chunk.sum(axis=(0, 1))
            label = "Order {}, Step {}".format(order, islice)
            if spectrum is not None:
                plt.figure(1)
                if xenergy is None:
                    plt.plot(spectrum, label=label)
                    plt.xlabel("MCA (channels)")
                else:
                    plt.plot(xenergy, spectrum, label=label)
                    plt.xlabel("MCA (keV)")
                plt.title("Exposure time = {} sec".format(time))
                plt.ylabel("Counts")
                if ylog:
                    plt.yscale("log")
                plt.legend()
            if img is not None:
                plt.figure()
                im = plt.imshow(chunk.sum(axis=-1), origin="lower", extent=extent)
                plt.title(label + " (looking upstream)")
                plt.xlabel("Synchrotron plane (cm)")
                plt.ylabel("Vertical direction (cm)")
                ax = plt.gca()
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(im, cax=cax)


class XrmcFile(object):
    def __init__(self):
        super(XrmcFile, self).__init__()

    @property
    def filename(self):
        if self.fileprefix and self.filesuffix:
            return self.fileprefix + "_" + self.filesuffix + ".dat"
        elif self.filesuffix:
            return self.filesuffix + ".dat"
        else:
            raise RuntimeError("Needs at least a file suffix")

    @property
    def filepath(self):
        return os.path.join(self.path, self.filename)

    @abstractproperty
    def body(self):
        pass

    @abstractproperty
    def header(self):
        pass

    @property
    def content(self):
        body = self.body
        if not body:
            return []
        lines = self.header
        if not lines:
            lines = []
        lines = ["; " + s for s in lines]
        lines.append("")
        lines += body
        lines += ["", "End"]
        return lines

    @property
    def empty(self):
        return not bool(self.body)

    def save(self, mode="w"):
        filepath = self.filepath
        logger.info("Saving " + filepath)
        localfs.Path(self.path).mkdir()
        localfs.Path(filepath).remove()
        lines = self.content
        if not lines:
            return
        with open(filepath, mode=mode) as f:
            for line in lines:
                f.write(line + "\n")

    @staticmethod
    def title(text):
        return ";%%%%%%%%%%%%%%%%%%%%% {} %%%%%%%%%%%%%%%%%%%%%".format(text.upper())


class LoopAxis(object):
    def __init__(self, parent, devices, suffix, steps=None, nsteps=1):
        """
        Args:
            parent(object):
            device(list(object)):
            suffix(str):
            steps(list(list(2-tuple))): ('Rotate', (...)), ('Translate', (...)), ('TranslateAll', (...)), ...
            nsteps(num): number of steps
        """
        self.transforms = []
        self.rtransforms = []
        self.nsteps = nsteps
        if steps is None:
            steps = [None] * len(devices)
        for device, step in zip(devices, steps):
            name = device.name + suffix
            transform = Transform(parent, name, device, step)
            self.transforms.append(transform)
            rtransform = transform.retour(nsteps)
            self.rtransforms.append(rtransform)

    def __repr__(self):
        return ";".join(repr(t) for t in self.transforms)

    def save(self, **kwargs):
        for transform in self.transforms:
            transform.save(**kwargs)
        for rtransform in self.rtransforms:
            rtransform.save(**kwargs)


class XrmcMain(XrmcFile):
    def __init__(self, path, fileprefix=None):
        self.path = path
        self.fileprefix = fileprefix
        self.devices = []
        self.loops = []  # first is the outer loop
        super(XrmcMain, self).__init__()

    def __repr__(self):
        s = ""
        lst = "\n ".join(map(repr, self.devices))
        if lst:
            s += "Devices:\n " + lst
        lst = "\n ".join([repr(ax) for ax in self.loops])
        if lst:
            if s:
                s += "\n"
            s += "Loops:\n " + lst
        return s

    @property
    def filesuffix(self):
        return "main"

    @property
    def header(self):
        return [
            "X-ray Monte Carlo (xrmc) input file",
            "Usage: xrmc {}".format(self.filename),
        ]

    @property
    def body(self):
        lines = []
        if self.devices:
            lines += ["", self.title("DEVICES USED FOR SIMULATION")]
            for device in self.devices:
                if not device.empty:
                    lines += ["", "Load " + device.filename]
        detectors = self.detectors()
        if detectors:
            lines += ["", self.title("START SIMULATIONS")]
            if self.loops:
                fmt = self.loop_suffix_format
                bfirst = True
                for idx, steps in self.loop_steps():
                    lines.append("")
                    if bfirst:
                        bfirst = False
                    else:
                        lines += [
                            "Load " + dstep.filename for step in steps for dstep in step
                        ]
                    for detector in detectors:
                        lines += [
                            "Run {}".format(detector.name),
                            "Save {} {} {}".format(
                                detector.name,
                                detector.img_type,
                                detector.reloutput(suffix=fmt.format(*idx)),
                            ),
                        ]
            else:
                for detector in detectors:
                    lines += [
                        "Run {}".format(detector.name),
                        "Save {} {} {}".format(
                            detector.name, detector.img_type, detector.reloutput()
                        ),
                    ]
        return lines

    def output_suffixes(self):
        if self.loops:
            fmt = self.loop_suffix_format
            return [fmt.format(*idx) for idx, _ in self.loop_steps()]
        else:
            return [""]

    def add_device(self, cls, *args, **kwargs):
        device = cls(self, *args, **kwargs)
        self.devices.append(device)
        return device

    def save(self, **kwargs):
        super(XrmcMain, self).save(**kwargs)
        for device in self.devices:
            device.save(**kwargs)
        for ax in self.loops:
            ax.save(**kwargs)

    def source(self, **kwargs):
        return self.find_device(cls=Source, **kwargs)

    def spectrum(self, **kwargs):
        return self.find_device(cls=Spectrum, **kwargs)

    def sample(self, **kwargs):
        return self.find_device(cls=Sample, **kwargs)

    def quadrics(self, **kwargs):
        return self.find_device(cls=Quadrics, **kwargs)

    def compositions(self, **kwargs):
        return self.find_device(cls=Compositions, **kwargs)

    def objects(self, **kwargs):
        return self.find_device(cls=Objects, **kwargs)

    def detectors(self, first=False, **kwargs):
        return self.find_device(cls=Detector, first=first, **kwargs)

    def find_device(self, cls=None, name=None, first=True):
        lst = []
        for device in self.devices:
            badd = True
            if cls is not None:
                badd &= isinstance(device, cls)
            if name is not None:
                badd &= device.name == name
            if badd:
                lst.append(device)
        if first:
            if lst:
                return lst[0]
            else:
                return None
        else:
            return lst

    def remove_device(self, cls=None, name=None):
        if cls is None and name is None:
            return
        lst = []
        for device in self.devices:
            bremove = True
            if cls is not None:
                bremove &= isinstance(device, cls)
            if name is not None:
                bremove &= device.name == name
            if not bremove:
                lst.append(device)
        self.devices = lst

    def simulate(self):
        if not installed():
            raise RuntimeError("xrmc is not installed")
        for detector in self.detectors():
            detector.removeoutput()
        localfs.Path(self.path)["output"].mkdir()
        _, err, returncode = execute(self.filename, cwd=self.path, stderr=True)
        success = returncode == 0
        if not success:
            if err:
                print("XRMC errors:")
                print(err)
            if "pulsetrain maximum reached" in err:
                self.create_empty_output()
                success = True
        return success

    def create_empty_output(self):
        ninteractions = self.sample().ninteractions
        for detector in self.detectors():
            detector.saveemptyoutput(ninteractions)

    def add_loop(self, devices, suffix, steps=None, nsteps=1):
        """
        First call specifies the outer loop and last call the inner loop

        Args:
            device(list(object)):
            suffix(str):
            steps(list(list(2-tuple))): ('Rotate', (...)), ('Translate', (...)), ('TranslateAll', (...)), ...
            nsteps(num): number of steps
        """
        self.loops.append(LoopAxis(self, devices, suffix, steps=steps, nsteps=nsteps))

    def removeloops(self):
        self.loops = []

    @property
    def loop_suffix_format(self):
        fmt = ""
        for ax in self.loops:
            fmt += "_{{:0{}d}}".format(max(len(str(ax.nsteps)), 1))
        return fmt

    def loop_steps(self):
        loops = self.loops
        loopidx = [list(range(ax.nsteps + 1)) for ax in loops]
        idx_prev = (-1,) * len(loops)
        for idx in itertools.product(*loopidx):
            idxsteps = []
            for ax, a, b in zip(loops, idx_prev, idx):
                if b > a:
                    idxsteps.append(ax.transforms)
                elif b < a:
                    idxsteps.append(ax.rtransforms)
            yield idx, idxsteps[::-1]
            idx_prev = idx


class XrmcChild(XrmcFile):
    def __init__(self, parent):
        self.parent = parent
        super(XrmcChild, self).__init__()

    @property
    def fileprefix(self):
        return self.parent.fileprefix

    @fileprefix.setter
    def fileprefix(self, value):
        self.parent.fileprefix = value

    @property
    def path(self):
        return self.parent.path

    @path.setter
    def path(self, value):
        self.parent.path = value

    @property
    def body(self):
        lines = []
        code = self.code
        if code:
            fmt = "{{:{}s}}".format(max(len(tpl[0]) for tpl in code if tpl))
            for tpl in code:
                if tpl:
                    var, comment, enabled = tpl
                    line = fmt.format(var)
                    if comment:
                        if enabled:
                            line += "  ; "
                        else:
                            line += " ; "
                        line += comment
                    if not enabled:
                        line = ";" + line
                else:
                    line = ""
                lines.append(line)
        return lines


class XrmcDevice(XrmcChild):

    TYPE = None

    def __init__(self, parent, name):
        super(XrmcDevice, self).__init__(parent)
        self.name = name

    def __repr__(self):
        return self.name

    @property
    def code(self):
        return [
            ("Newdevice " + self.TYPE, "Object type", True),
            (self.name, "Object name", True),
        ]

    @property
    def filesuffix(self):
        return self.name

    @filesuffix.setter
    def filesuffix(self, value):
        self.name = value.lower()


class XrmcReferenceFrame(object):
    """
    The XRMC sample reference frame:

    * Y-axis: beam direction
    * X-axis: horizontal and pointing to the left when looking upstream
    * Z-axis: vertical and pointing upwards

    The User sample reference frame:

    * Z-axis: beam direction
    * X-axis: horizontal and pointing to the right when looking upstream
    * Y-axis: vertical and pointing upwards

    These reference frames are related by x', y', z' = -x, z, y.

    Additionally, the rotation angle in the User and XRMC reference frame
    is in the opposite direction.
    """

    def __init__(self):
        super(XrmcReferenceFrame, self).__init__()

    @staticmethod
    def xyzUserToXrmc(x, y, z):
        return -x, z, y

    def rotationUserToXrmc(self, x0, y0, z0, ux, uy, uz, angle):
        return (
            self.xyzUserToXrmc(x0, y0, z0) + self.xyzUserToXrmc(ux, uy, uz) + (-angle,)
        )

    @staticmethod
    def xrmcArrayInput(name, comment, a, identity=False):
        return " ".join([name] + list(map(str, a))), comment, not identity

    def xrmcTranslationInput(self, name, item, *args, **kw):
        """Defined in USER frame, return in XRMC frame
        """
        arr = self.xyzUserToXrmc(*args, **kw)
        comment = "{} in the sample frame (x, y, z; cm)".format(item)
        return self.xrmcArrayInput(
            name, comment, arr, identity=all(x == 0 for x in arr)
        )

    def xrmcRotationInput(self, name, *args, **kw):
        """Defined in USER frame, return in XRMC frame
        """
        arr = self.rotationUserToXrmc(*args, **kw)
        comment = "rotation around axis passing through (x, y, z; cm), with direction (u,v,w; cm) and rotation angle (deg)"
        return self.xrmcArrayInput(name, comment, arr, identity=arr[-1] == 0)

    def xrmcTranslationMatrix(self, *args, **kw):
        """Defined in USER frame, return in XRMC frame
        """
        M = np.eye(4)
        M[0:3, 3] = self.xyzUserToXrmc(*args, **kw)
        return M

    def xrmcRotationMatrix(self, *args, **kw):
        """Defined in USER frame, return in XRMC frame
        """
        arr = self.rotationUserToXrmc(*args, **kw)
        angle = np.radians(arr[-1])
        if angle == 0:
            return np.eye(4)
        arr = np.array(arr)
        pt = arr[:3]
        vec = arr[3:6]
        rotvec = vec * (angle / sum(vec * vec))
        M1 = np.eye(4)
        M1[:3, 3] = -pt
        M2 = np.eye(4)
        M2[:3, :3] = Rotation.from_rotvec(rotvec).as_matrix()
        M3 = np.eye(4)
        M3[:3, 3] = pt
        return M3.dot(M2.dot(M1))

    @property
    def userToXrmc(self):
        return np.array([[-1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])

    @property
    def xrmcToUser(self):
        return np.array([[-1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])

    def applyTransform(self, coord):
        """
        Args:
            coord(array): user coordinates
        Returns:
            array: xrmc coordinates
        """
        coord = np.atleast_2d(coord)
        coord = np.vstack([coord, np.ones(coord.shape[1])])
        coord = self.userToXrmc.dot(coord.T).T
        return coord[:, :-1] / coord[:, -1]


class XrmcPositionalDevice(XrmcDevice, XrmcReferenceFrame):
    """
    Device reference point is positioned at spherical coordinates [r, polar, azimuth].
    Device surface goes through the reference point.
    Device surface-center coincides with reference point (unless offsets are not zero).
    """

    def __init__(
        self,
        parent,
        name,
        radius=0,
        polar=0,
        azimuth=0,
        hoffset=0,
        voffset=0,
        fixed=True,
    ):
        super(XrmcPositionalDevice, self).__init__(parent, name)
        self.radius = radius
        self.polar = polar
        self.azimuth = azimuth
        self.hoffset = hoffset
        self.voffset = voffset
        self.fixed = fixed

    @property
    def userToXrmc(self):
        """The user frame of the device (i.e. before moving the device)
        """
        M0 = super().userToXrmc
        M1 = self.xrmcTranslationMatrix(self.hoffset, self.voffset, self.radius)
        M2 = self.xrmcRotationMatrix(0, 0, 0, 0, 1, 0, self.polar)
        M3 = self.xrmcRotationMatrix(0, 0, 0, 0, 0, 1, self.azimuth)
        return M3.dot(M2.dot(M1).dot(M0))

    @property
    def xrmcToUser(self):
        M0i = super().xrmcToUser
        M1i = self.xrmcTranslationMatrix(-self.hoffset, -self.voffset, -self.radius)
        M2i = self.xrmcRotationMatrix(0, 0, 0, 0, 1, 0, -self.polar)
        M3i = self.xrmcRotationMatrix(0, 0, 0, 0, 0, 1, -self.azimuth)
        return M0i.dot(M1i.dot(M2i).dot(M3i))

    @property
    def quadricChangeOfFrame(self):
        """Quadric already in XRMC frame
        """
        M1i = self.xrmcTranslationMatrix(-self.hoffset, -self.voffset, -self.radius)
        M2i = self.xrmcRotationMatrix(0, 0, 0, 0, 1, 0, self.polar)
        M3i = self.xrmcRotationMatrix(0, 0, 0, 0, 0, 1, self.azimuth)
        return M1i.dot(M2i.dot(M3i))

    def xmrcTransformedQuadric(self, name, typ, params):
        """Quadric params already in XRMC frame
        """
        idx = ([0, 1, 2, 3, 1, 2, 3, 2, 3, 3], [0, 0, 0, 0, 1, 1, 1, 2, 2, 3])
        if typ == "Plane":
            A = quadrics.plane(params[:3], params[3:])
        elif typ == "Ellipsoid":
            A = quadrics.ellipsoid(params[:3], params[3:])
        elif typ.startswith("Cylinder"):
            axis = ["X", "Y", "Z"].index(typ[-1])
            A = quadrics.cylinder(params[:3], params[3:], axis)
        else:
            A = np.zeros((4, 4))
            A[idx] = params
        A = quadrics.change_of_frame(A, self.quadricChangeOfFrame)
        return A[idx].tolist()

    def xrmcTransformInput(self, item, all=True):
        """
        Transformation of point [0, 0, 0] (cartesian) to [r, polar, azimuth] (spherical)
        """
        # Spherical: [r, polar, azimuth]
        # Cartesian: RotZ(azimuth).RotY(polar).Trn([0, 0, r]).[0, 0, 0]
        if all:
            rname = "RotateAll"
            tname = "TranslateAll"
        else:
            rname = "Rotate"
            tname = "Translate"
        return [
            self.xrmcTranslationInput(
                tname, item, self.hoffset, self.voffset, self.radius
            ),
            self.xrmcRotationInput(rname, 0, 0, 0, 0, 1, 0, self.polar),
            self.xrmcRotationInput(rname, 0, 0, 0, 0, 0, 1, self.azimuth),
        ]

    def xrmcRotateInput(self, all=True):
        """
        Transformation of point [0, 0, r] (cartesian) to [r, polar, azimuth] (spherical)
        """
        # Spherical: [r, polar, azimuth]
        # Cartesian: RotZ(azimuth).RotY(polar).[0, 0, r]
        if all:
            name = "RotateAll"
        else:
            name = "Rotate"
        return [
            self.xrmcRotationInput(name, 0, 0, 0, 0, 1, 0, self.polar),
            self.xrmcRotationInput(name, 0, 0, 0, 0, 0, 1, self.azimuth),
        ]

    def xrmcPositionInput(self, item):
        """
        Position device at [r, polar, azimuth] (spherical)
        """
        return self.xrmcTranslationInput(
            "X", item, self.hoffset, self.voffset, self.radius
        )

    def add_rotationloop(self, suffix, x, y, z, u, v, w, step, nsteps):
        """
        Rotate around axis going through [x, y, z] with direction [u, v, w]
        """
        params = x, y, z, u, v, w, step
        if self.fixed:
            raise RuntimeError("{} is fixed".format(repr(self)))
        devices = []
        steps = []
        # Rotate the device quadrics (e.g. sample, detector filters)
        devices.append(self.quadrics)
        steps.append([("RotateAll", params)])
        if not isinstance(self, Sample):
            # When rotating a device, the sample should be fixed
            sample = self.parent.sample()
            if not sample.fixed:
                raise RuntimeError("{} must be fixed".format(repr(sample)))
            # Rotate the device itself with the quadrics
            devices.append(self)
            steps.append([("Rotate", params)])
        self.parent.add_loop(devices, suffix, steps=steps, nsteps=nsteps)

    def add_Xrotationloop(self, step, nsteps):
        self.add_rotationloop("_Xrot", 0, 0, 0, 1, 0, 0, step, nsteps)

    def add_Yrotationloop(self, step, nsteps):
        self.add_rotationloop("_Yrot", 0, 0, 0, 0, 1, 0, step, nsteps)

    def add_Zrotationloop(self, step, nsteps):
        self.add_rotationloop("_Zrot", 0, 0, 0, 0, 0, 1, step, nsteps)

    def add_polarrotationloop(self, step, nsteps):
        """
        Keeps azimuth constant modulo 180 degrees
        """
        angle = np.radians(self.azimuth + 90)
        u = np.cos(angle)
        v = np.sin(angle)
        self.add_rotationloop("_polarrot", 0, 0, 0, u, v, 0, step, nsteps)

    def add_azimuthalrotationloop(self, step, nsteps):
        """
        Keeps polar angle constant
        """
        self.add_rotationloop("_azimuthrot", 0, 0, 0, 0, 0, 1, step, nsteps)

    def add_layer(
        self,
        material=None,
        thickness=None,
        dhor=None,
        dvert=None,
        ohor=0,
        overt=0,
        **kwargs
    ):
        if not material:
            material = "Vacuum"
        ilayer, names = self.quadrics.add_layer(
            thickness, dhor, dvert, ohor=ohor, overt=overt, pdevice=self, **kwargs
        )
        matname = self.compositions.add_material(material)
        objname = "{}_layer{}".format(self.name, ilayer)
        self.objects.add(objname, material, names)
        return objname, matname

    @property
    def quadrics(self):
        return self.parent.quadrics()

    @property
    def objects(self):
        return self.parent.objects()

    @property
    def compositions(self):
        return self.parent.compositions()


class Transform(XrmcChild, XrmcReferenceFrame):
    def __init__(self, parent, name, device, transformations):
        super(Transform, self).__init__(parent)
        self.filesuffix = name
        self.device = device
        self.transformations = transformations
        for cmd, parameters in transformations:
            if cmd.startswith("Rotate"):
                if len(parameters) != 7:
                    raise ValueError("A rotation needs 7 parameters: x,y,z,u,v,w,angle")
            elif cmd.startswith("Translate"):
                if len(parameters) != 3:
                    raise ValueError("A translation needs 3 parameters: dx,dy,dz")
            else:
                raise ValueError("Unknown transformation {}".format(repr(cmd)))

    def __repr__(self):
        return ",".join(cmd for cmd, _ in self.transformations) + "({})".format(
            self.device
        )

    @property
    def header(self):
        return [
            "Transformation of device {}.".format(repr(self.device.name)),
            "Executed every time this file is loaded.",
        ]

    @property
    def code(self):
        lines = [("Device {}".format(self.device.name), "device name", True)]
        for cmd, parameters in self.transformations:
            if cmd.startswith("Rotate"):
                lines.append(self.xrmcRotationInput(cmd, *parameters))
            elif cmd.startswith("Translate"):
                lines.append(self.xrmcTranslationInput(cmd, "", *parameters))
            else:
                raise NotImplementedError(cmd)
        return lines

    def retour(self, nsteps):
        rtransformations = []
        for cmd, parameters in self.transformations:
            if cmd.startswith("Rotate"):
                parameters = tuple(parameters[:-1]) + (-nsteps * parameters[-1],)
            elif cmd.startswith("Translate"):
                parameters = (
                    -nsteps * parameters[0],
                    -nsteps * parameters[1],
                    -nsteps * parameters[2],
                )
            else:
                raise NotImplementedError(cmd)
            rtransformations.append((cmd, parameters))
        return self.__class__(
            self.parent, self.filesuffix + "_retour", self.device, rtransformations
        )


class Source(XrmcPositionalDevice):

    TYPE = "source"

    def __init__(self, parent, name, distance=None, beamsize=None):
        super(Source, self).__init__(
            parent, name, radius=distance, polar=180, azimuth=0, hoffset=0, voffset=0
        )
        self.beamsize = beamsize

    @property
    def distance(self):
        return self.radius

    @distance.setter
    def distance(self, value):
        self.radius = value

    @property
    def position(self):
        return -self.distance

    @property
    def header(self):
        return [
            "Source position/orientation file",
            "",
            "The sample reference frame is defined so that the beam travels along",
            "the Y-axis, the Z-axis points upwards and the X-axis is parallel with",
            "the synchrotron plane, pointing to the left when looking upstream.",
            "",
            "In the local reference from of the source, the beam",
            "travels along the Z-axis.",
            "",
            "The source is positioned at y=-distance and has a divergence of",
            "(dx, dy) and size (dx, dy, dz) in the local reference frame.",
        ]

    @property
    def divergence(self):
        return np.arctan2(self.beamsize * 0.5, self.distance)

    @property
    def code(self):
        lines = super(Source, self).code
        div = self.divergence
        position = self.xrmcTranslationInput("X", "source origin", 0, 0, -self.distance)

        # Source reference frame:
        #   XY-plane: Electric field vector [Ex, Ey]
        #   Z: beam direction
        # X-axis is user frame: horizontal, pointing to the right when looking upstream
        # Z-axis is user frame: horizontal, pointing downstream
        xdirection = self.xrmcTranslationInput("ui", "source X-axis", 1, 0, 0)
        zdirection = self.xrmcTranslationInput("uk", "source Z-axis", 0, 0, 1)
        lines += [
            (
                "SpectrumName {}".format(self.parent.spectrum().name),
                "Spectrum input device name",
                True,
            ),
            position,
            xdirection,
            zdirection,
            (
                "Divergence {} {}".format(div, div),
                "beam divergency (dx, dy; rad)",
                True,
            ),
            ("Size 0.0 0.0 0.0", "source size (sx, sy, sz; cm)", True),
        ]
        return lines


class Spectrum(XrmcDevice):

    TYPE = "spectrum"

    def __init__(
        self,
        parent,
        name,
        lines=None,
        polarization=True,
        polarization_angle=0,
        multiplicity=1,
    ):
        super(Spectrum, self).__init__(parent, name)
        self.polarization = polarization
        self.lines = lines
        self.multiplicity = multiplicity
        self.polarization_angle = polarization_angle

    @property
    def header(self):
        return [
            "Spectrum file",
            "",
            "Discrete X-ray source spectrum + continuum",
            "",
            "Polarization is either unpolarized or linear polarized.",
            "The orientation of the electric field vector is defined",
            "by (Ex, Ey) in the source reference frame (see source).",
        ]

    @property
    def code(self):
        lines = super(Spectrum, self).code
        m = int(np.round(self.multiplicity))
        lines += [
            (
                "PolarizedFlag {:d}".format(self.polarization),
                "unpolarized/polarized beam (0/1)",
                True,
            ),
            ("LoopFlag 1", "0: extract random energies on the whole spectrum", True),
            ("", "1: loop on all lines and sampling points", True),
            (
                "ContinuousPhotonNum {:d}".format(m),
                "Multiplicity of events for each spectrum interval",
                True,
            ),
            (
                "LinePhotonNum {:d}".format(m),
                "Multiplicity of events for each spectrum line",
                True,
            ),
            ("RandomEneFlag 1", "enable random energy on each interval (0/1)", True),
            ("Lines", "discrete energy lines of the spectrum", True),
            (str(len(self.lines)), "Number of lines in the spectrum", True),
        ]
        cos = np.cos(self.polarization_angle)
        sin = np.sin(self.polarization_angle)
        for energy, sigma, mult in self.lines:
            lines.append(
                (
                    "{} {} {} {}".format(energy, sigma, mult * cos, mult * sin),
                    "Energy (keV), sigma (keV), Ex (ph/s), Ey (ph/s)",
                    True,
                )
            )
        return lines

    @property
    def emax(self):
        return max(energy for energy, sigma, intensity in self.lines)

    @property
    def totalflux(self):
        return sum(intensity for energy, sigma, intensity in self.lines)


class Detector(XrmcPositionalDevice):

    TYPE = "detectorarray"

    def __init__(
        self,
        parent,
        name,
        distance=None,
        poissonnoise=True,
        ebinsize=None,
        forcedetect=False,
        multiplicity=None,
        emin=None,
        emax=None,
        time=1,
        pixelshape="elliptical",
        as_lines=True,
        **kwargs
    ):
        super(Detector, self).__init__(parent, name, radius=distance, **kwargs)
        self.emin = emin
        self.emax = emax
        self.ebinsize = ebinsize
        self.poissonnoise = poissonnoise
        self.forcedetect = forcedetect
        self.multiplicity = multiplicity
        self.time = time
        self.pixelshape = pixelshape.lower()
        self.as_lines = as_lines

    @property
    def distance(self):
        return self.radius

    @distance.setter
    def distance(self, value):
        self.radius = value

    @property
    def position(self):
        return self.distance

    @property
    def header(self):
        return [
            "Detector array parameter file",
            "",
            "The detector reference frame (before transformation) is defined",
            "so that the Z-axis points upstream, the X-axis (detector image rows)",
            "points upwards and the Y-axis (detector image columns) points to the",
            "right when looking upstream.",
            "",
            "A detector can have 1 or more elliptical or square pixels.",
        ]

    @property
    def code(self):
        lines = super(Detector, self).code
        energydispersive = self.energydispersive
        m = int(np.round(self.multiplicity))
        if energydispersive:
            pixeltype = 2
        else:
            pixeltype = 0

        position = self.xrmcPositionInput("detector origin")
        rotations = self.xrmcRotateInput(all=False)
        # Detector reference frame:
        #   X-axis: rows
        #   Y-axis: columns
        #   Z-axis: surface normal (should point to the sample)
        # X-axis is user frame: vertical, pointing upward
        # Z-axis is user frame: horizontal, pointing upstream
        xdirection = self.xrmcTranslationInput("ui", "detector X-axis (rows)", 0, 1, 0)
        zdirection = self.xrmcTranslationInput(
            "uk", "detector Z-axis (surface normal)", 0, 0, -1
        )
        lines += [
            (
                "SourceName {}".format(self.parent.sample().name),
                "Source input device name",
                True,
            ),
            (
                "NPixels {} {}".format(*self.dims[::-1]),
                "Pixel number (nx x ny, ncol x nrow)",
                True,
            ),
            (
                "PixelSize {} {}".format(*self.pixelsize[::-1]),
                "Pixel Size (sx x sy, scol x srow; cm)",
                True,
            ),
            (
                "Shape {:d}".format(self.elliptical),
                "Pixel shape (0 rectangular, 1 elliptical)",
                True,
            ),
            (
                "dOmegaLim {}".format(2 * np.pi),
                "Limit solid angle (default = 2*PI)",
                False,
            ),
            position,
            xdirection,
            zdirection,
            ("ExpTime {}".format(self.time), "Exposure time (sec)", True),
            (
                "PhotonNum {:d}".format(m),
                "Multiplicity of simulated events per pixel",
                True,
            ),
            ("RandomPixelFlag 1", "Enable random point on pixels (0/1)", True),
            (
                "PoissonFlag {:d}".format(self.poissonnoise),
                "Enable Poisson statistic on pix. counts (0/1)",
                True,
            ),
            ("RoundFlag 0", "Round pixel counts to integer (0/1)", True),
            ("AsciiFlag 0", "binary(0) or ascii(1) file format", False),
            (
                "ForceDetectFlag {:d}".format(self.forcedetect),
                "MC variance reduction (0/1)",
                True,
            ),
            ("HeaderFlag 1", "Use header in output file (0/1)", True),
            ("PixelType {:d}".format(pixeltype), "Pixel content type:", True),
            ("", "0: fluence,      1: energy fluence,", True),
            ("", "2: fluence(E),   3: energy fluence(E)", True),
            ("Emin {}".format(self.emin), "Emin", energydispersive),
            ("Emax {}".format(self.emax), "Emax", energydispersive),
            ("NBins {}".format(self.nbins), "NBins", energydispersive),
            (
                "SaturateEmin 0",
                "Saturate energies lower than Emin (0/1)",
                energydispersive,
            ),
            (
                "SaturateEmax 0",
                "Saturate energies greater than Emax (0/1)",
                energydispersive,
            ),
        ]
        lines += rotations
        return lines

    def reloutput(self, suffix=None):
        if not suffix:
            suffix = ""
        return "output/{}{}.dat".format(self.name, suffix)

    def absoutput(self, suffix=None):
        if not suffix:
            suffix = ""
        return os.path.join(self.outpath, "{}{}.dat".format(self.name, suffix))

    @property
    def outpath(self):
        return os.path.join(self.path, "output")

    @property
    def energydispersive(self):
        return bool(self.nbins)

    @property
    def elliptical(self):
        return self.pixelshape == "elliptical"

    @property
    def img_type(self):
        if self.as_lines:
            return "Image"
        else:
            return "ConvolutedImage"

    @property
    def emax(self):
        if self._emax is None:
            return self.parent.spectrum().emax + 1
        else:
            return self._emax

    @emax.setter
    def emax(self, value):
        self._emax = value

    @property
    def emin(self):
        if self._emin is None:
            return 0
        else:
            return self._emin

    @emin.setter
    def emin(self, value):
        self._emin = value

    @property
    def ebinsize(self):
        if self.energydispersive:
            return (self.emax - self.emin) / (self.nbins - 1.0)
        else:
            return 0

    @ebinsize.setter
    def ebinsize(self, value):
        if bool(value):
            self.nbins = max(
                int(np.round((self.emax - self.emin) / float(value) + 1)), 1
            )
        else:
            self.nbins = 1

    @property
    def mcagain(self):
        return self.ebinsize

    @property
    def mcazero(self):
        return 0

    def result(self):
        if self.parent.loops:
            suffix = "_*"
        else:
            suffix = None
        filename = self.absoutput(suffix=suffix)
        dirname = os.path.dirname(filename)
        basename, ext = os.path.splitext(os.path.basename(filename))
        return loadxrmcresult(dirname, basename, ext)

    def removeoutput(self):
        filename = self.absoutput()
        if os.path.isfile(filename):
            os.remove(filename)

    def saveemptyoutput(self, ninteractions):
        energydispersive = self.energydispersive
        if energydispersive:
            pixeltype = 2
        else:
            pixeltype = 0
        nrow, ncol = self.dims
        srow, scol = self.pixelsize
        header = np.array(
            (
                ninteractions,
                ncol,
                nrow,
                scol,
                srow,
                self.time,
                pixeltype,
                self.nbins,
                self.emin,
                self.emax,
            ),
            dtype=headerdt,
        )
        suffixes = self.parent.output_suffixes()
        for suffix in suffixes:
            saveemptyresult(self.absoutput(suffix=suffix), header)

    def add_layer(self, thickness=None, **kwargs):
        """Add layer in front of the detector
        """
        thickness = -abs(thickness)
        super(Detector, self).add_layer(thickness=thickness, **kwargs)


class AreaDetector(Detector):

    TYPE = "detectorarray"

    def __init__(
        self,
        parent,
        name,
        pixelsize=None,
        dims=None,
        hpoffset=None,
        vpoffset=None,
        **kwargs
    ):
        self.dims = dims
        self.pixelsize = pixelsize
        super(AreaDetector, self).__init__(parent, name, **kwargs)
        if hpoffset is not None:
            self.hpoffset = hpoffset
        if vpoffset is not None:
            self.vpoffset = vpoffset

    @property
    def hpoffset(self):
        return self.hoffset / float(self.pixelsize[0])

    @hpoffset.setter
    def hpoffset(self, value):
        self.hoffset = value * self.pixelsize[0]

    @property
    def vpoffset(self):
        return self.voffset / float(self.pixelsize[1])

    @vpoffset.setter
    def vpoffset(self, value):
        self.voffset = value * self.pixelsize[1]

    @property
    def size(self):
        return self.dims[0] * self.pixelsize[0], self.dims[1] * self.pixelsize[1]


class SingleElementDetector(Detector):

    TYPE = "detectorarray"

    def __init__(self, parent, name, activearea=None, **kwargs):
        self.dims = 1, 1
        self.activearea = activearea
        super(SingleElementDetector, self).__init__(parent, name, **kwargs)

    @property
    def pixelsize(self):
        if self.elliptical:
            a = 2 * (self.activearea / np.pi) ** 0.5
        else:
            a = self.activearea ** 0.5
        return a, a


class SDD(SingleElementDetector):

    TYPE = "detectorconvolute"

    def __init__(
        self,
        parent,
        name,
        material=None,
        thickness=None,
        windowmaterial=None,
        windowthickness=0,
        dtfrac=None,
        countrate=None,
        pulseproctime=0,
        noise=None,
        fano=None,
        as_lines=False,
        **kwargs
    ):
        self.noise = noise
        self.fano = fano
        self.set_pulseproctime(dtfrac, countrate, pulseproctime=pulseproctime)
        super(SDD, self).__init__(parent, name, as_lines=as_lines, **kwargs)
        material = self.parent.compositions().add_material(material)
        windowmaterial = self.parent.compositions().add_material(windowmaterial)
        self.material = material
        self.thickness = thickness
        self.windowmaterial = windowmaterial
        self.windowthickness = windowthickness

    @property
    def code(self):
        lines = super(SDD, self).code
        lines += [
            (
                "CompositionName {}".format(self.parent.compositions().name),
                "Composition input device name",
                True,
            ),
            (
                "CrystalPhase {}".format(self.material),
                "Detector crystal material",
                True,
            ),
            (
                "CrystalThickness {}".format(self.thickness),
                "Detector crystal thickness",
                True,
            ),
            ("WindowPhase {}".format(self.windowmaterial), "Window material", True),
            (
                "WindowThickness {}".format(self.windowthickness),
                "Window thickness",
                True,
            ),
            ("Noise {}".format(self.noise), "Detector energy noise (keV)", True),
            (
                "FanoFactor {}".format(self.fano),
                "Detector fano noise (dimensionless)",
                True,
            ),
            (
                "PulseWidth {}".format(self.pulseproctime),
                "Time to process one pulse (sec)",
                bool(self.pulseproctime),
            ),
        ]
        return lines

    def set_pulseproctime(self, dtfrac, countrate, pulseproctime=0):
        # Assume extendable deadtime
        if dtfrac and countrate:
            self.pulseproctime = -np.log(1 - dtfrac) / countrate
        else:
            self.pulseproctime = pulseproctime


class Sample(XrmcPositionalDevice):

    TYPE = "sample"

    def __init__(self, parent, name, multiplicity=(1, 1), fixed=False, **kwargs):
        super(Sample, self).__init__(parent, name, fixed=fixed, **kwargs)
        self.multiplicity = multiplicity

    @property
    def ninteractions(self):
        return len(self.multiplicity)

    @property
    def header(self):
        return [
            "Sample parameters file",
            "",
            "maximum scattering order",
            "    0: transmission",
            "    1: first order scattering or fluorescence emission",
            "    2: second order scattering or fluorescence emission",
        ]

    @property
    def code(self):
        lines = super(Sample, self).code
        lines += [
            (
                "SourceName {}".format(self.parent.source().name),
                "Source input device name",
                True,
            ),
            (
                "Geom3DName {}".format(self.parent.objects().name),
                "Geom3d input device name",
                True,
            ),
            (
                "CompName {}".format(self.parent.compositions().name),
                "Composition input device name",
                True,
            ),
            ("WeightedStepLength 1", "Weighted step length (0/1)", True),
            (
                "FluorFlag {:d}".format(len(self.multiplicity) > 1),
                "Activate Fluorescence (0/1)",
                True,
            ),
            (
                "ScattOrderNum {:d}".format(len(self.multiplicity) - 1),
                "Maximum scattering order",
                True,
            ),
        ]
        for i, mi in enumerate(self.multiplicity):
            m = str(int(np.round(mi)))
            lines.append(
                (m, "Multiplicity of simulated events for order {}".format(i), True)
            )
        return lines


class Objects(XrmcDevice):

    TYPE = "geom3d"

    def __init__(self, parent, name):
        super(Objects, self).__init__(parent, name)
        self.clear()

    @property
    def header(self):
        return [
            "3D object geometric description file",
            "",
            "Makes objects with quadrics and compositions.",
        ]

    @property
    def code(self):
        if not self._dict:
            return []
        lines = super(Objects, self).code
        lines += [
            (
                "QArrName {}".format(self.parent.quadrics().name),
                "Quadric array input device name",
                True,
            ),
            (
                "CompName {}".format(self.parent.compositions().name),
                "Composition input device name",
                True,
            ),
        ]
        for name, (compoundin, compoundout, quadricnames) in self._dict.items():
            compoundin = materialname(compoundin)
            compoundout = materialname(compoundout)
            lines += [
                None,
                ("Object {}".format(name), "", True),
                ("{} {}".format(compoundin, compoundout), "Phase in, phase out", True),
                (str(len(quadricnames)), "Number of quadrics", True),
                (" ".join(quadricnames), "Quadric names", True),
            ]
        return lines

    def add(self, name, compoundin, quadricnames):
        try:
            compoundout = self.parent.compositions().atmosphere
        except AttributeError:
            compoundout = None
        self._dict[name] = [compoundin, compoundout, quadricnames]

    def change_atmosphere(self, compoundout):
        for lst in self._dict.values():
            lst[1] = compoundout

    def clear(self):
        self._dict = {}


class QuadricsGroup(MutableMapping):
    def __init__(self, pdevice, **kw):
        self.pdevice = pdevice
        self._qdict = OrderedDict(**kw)

    def __str__(self):
        return str(self._qdict)

    def __repr__(self):
        return repr(self._qdict)

    def __getitem__(self, key):
        return self._qdict[key]

    def __setitem__(self, key, value):
        self._qdict[key] = value

    def __delitem__(self, key):
        del self._qdict[key]

    def __len__(self):
        return len(self._qdict)

    def __iter__(self):
        return iter(self._qdict)


class Quadrics(XrmcDevice):

    TYPE = "quadricarray"

    def __init__(self, parent, name, **kwargs):
        """
        Args:
            parent
            name
        """
        super(Quadrics, self).__init__(parent, name, **kwargs)
        self.clear()

    def clear(self):
        self._qdict = OrderedDict()

    def group(self, pdevice):
        """Add group of quadrics associated to a positional device (e.g. sample, detector)

        Args:
            pdevice(XmrcPositionalDevice)
        """
        name = pdevice.name
        grp = self._qdict.get(name)
        if grp is None:
            self._qdict[name] = grp = QuadricsGroup(pdevice)
        return grp

    @property
    def nonfixed_group(self):
        for o in self._qdict.values():
            if not o.pdevice.fixed:
                return o

    @property
    def header(self):
        return [
            "Quadric array file",
            "",
            "A quadric (a.k.a. quadratic surface) is a parametric surface (plane, cylinder, sphere, ...):",
            "  x.A.x^T = 0   (where x homogeneous coordinates)",
            "For internal space x.A.x^T < 0 and for external space x.A.x^T > 0",
        ]

    @property
    def code(self):
        if not self._qdict:
            return []
        lines = super(Quadrics, self).code
        if sum(not o.pdevice.fixed for o in self._qdict.values()) > 1:
            raise RuntimeError("Only one non-fixed device with quadrics allowed")
        for grp in self._qdict.values():
            pdevice = grp.pdevice
            for name, typ, params in grp.values():
                if pdevice.fixed:
                    params = pdevice.xmrcTransformedQuadric(name, typ, params)
                    name += " BlockTransformAll"
                    comment = (
                        typ
                        + " A0*x^2 + 2*A1*x*y + 2*A2*x*z + 2*A3*x + A4*y^2 + 2*A5*y*z + 2*A6*y + A7*z^2 + 2*A8*z + A9 = 0"
                    )
                    typ = "Quadric"
                else:
                    if typ == "Plane":
                        comment = (
                            "Plane through (x0,y0,z0) with normal vector (nx,ny,nz)"
                        )
                    elif typ.startswith("Ellipsoid"):
                        axis = typ[-1]
                        comment = "Ellipsoid with center (xa, xb, xc) in the perpendicular plane. (Ra, Rb, Rc) are the elliptical semi-axes."
                    elif typ.startswith("Cylinder"):
                        axis = typ[-1]
                        comment = "Cylinder with main axis along {} and with coordinates (xa, xb) in the perpendicular plane. (Ra, Rb) are the elliptical semi-axes.".format(
                            axis
                        )
                    elif typ.startswith("Quadric"):
                        comment = "Equation A0*x^2 + 2*A1*x*y + 2*A2*x*z + 2*A3*x + A4*y^2 + 2*A5*y*z + 2*A6*y + A7*z^2 + 2*A8*z + A9 = 0"
                lines.append(("", "", True))
                lines.append(("{} {}".format(typ, name), "", True))
                lines.append((" ".join(list(map(str, params))), comment, True))
        lines.append(("", "", True))
        grp = self.nonfixed_group
        if grp is not None:
            lines += grp.pdevice.xrmcTransformInput(repr(grp.pdevice.name))
        return lines

    def add_plane(self, objname, x0, y0, z0, nx, ny, nz, pdevice=None):
        """Already in XMRC reference frame
        """
        grp = self.group(pdevice)
        if objname in grp:
            raise ValueError("Quadric {} already exists".format(repr(objname)))
        params = [x0, y0, z0, nx, ny, nz]
        name = grp.pdevice.name + "_" + objname
        grp[objname] = (name, "Plane", params)
        return name

    def add_cylinder(self, objname, axis, xa, xb, Ra, Rb, pdevice=None):
        """Already in XMRC reference frame
        """
        grp = self.group(pdevice)
        if objname in grp:
            raise ValueError("Quadric {} already exists".format(repr(objname)))
        axis = axis.upper()
        if axis not in ["X", "Y", "Z"]:
            raise ValueError("Axis must be 'X', 'Y' or 'Z'")
        params = [xa, xb, Ra, Rb]
        name = grp.pdevice.name + "_" + objname
        grp[objname] = (name, "Cylinder" + axis, params)
        return name

    def add_ellipsoid(self, objname, xa, xb, xc, Ra, Rb, Rc, pdevice=None):
        """Already in XMRC reference frame
        """
        grp = self.group(pdevice)
        if objname in grp:
            raise ValueError("Quadric {} already exists".format(repr(objname)))
        params = [xa, xb, xc, Ra, Rb, Rc]
        name = grp.pdevice.name + "_" + objname
        grp[objname] = (name, "Ellipsoid", params)
        return name

    def add_quadric(
        self, objname, A11, A12, A13, A14, A22, A23, A24, A33, A34, A44, pdevice=None
    ):
        """Already in XMRC reference frame
        """
        grp = self.group(pdevice)
        if objname in grp:
            raise ValueError("Quadric {} already exists".format(repr(objname)))
        params = [A11, A12, A13, A14, A22, A23, A24, A33, A34, A44]
        name = grp.pdevice.name + "_" + objname
        grp[objname] = (name, "Quadric", params)
        return name

    def add_layer(
        self,
        thickness,
        dhor,
        dvert,
        ohor=0,
        overt=0,
        surface=0,
        spacing=1e-10,
        pdevice=None,
    ):
        """
        Layer perpendicular to the beam direction (before transformation).
        Beam runs through the geometric center when `ohor == overt == 0`.
        Layer starts at the previous `layer + spacing` and ends at
        `layer + spacing + thickness`.

        Args:
            thickness: can be negative
            dhor: length along the horizontal axis
            dvert: length along the vertical axis
            ohor: shift along horizontal axis
            overt: shift along vertical axis
            surface: position of the surface on the perpendicular axis
            spacing: layers cannot be touching
        """
        ilayer = -1
        if thickness * spacing < 0:
            spacing = -spacing
        for name, clsname, params in self.group(pdevice).values():
            if clsname != "Plane":
                continue
            m = re.search(r"layer(\d+)b", name)
            if m:
                i = int(m.groups()[0])
                if i > ilayer:
                    ilayer = i
                    surface = params[1]  # end of previous layer
        ilayer += 1
        ta = surface + spacing
        tb = ta + thickness
        if thickness < 0:
            ta, tb = tb, ta
        names = []
        names.append(
            self.add_plane(
                "layer{}a".format(ilayer), 0, ta, 0, 0, -1, 0, pdevice=pdevice
            )
        )
        names.append(
            self.add_plane(
                "layer{}b".format(ilayer), 0, tb, 0, 0, 1, 0, pdevice=pdevice
            )
        )
        names.append(
            self.add_plane(
                "blayer{}_xmin".format(ilayer),
                ohor - dhor / 2,
                0,
                0,
                -1,
                0,
                0,
                pdevice=pdevice,
            )
        )
        names.append(
            self.add_plane(
                "blayer{}_xmax".format(ilayer),
                ohor + dhor / 2,
                0,
                0,
                1,
                0,
                0,
                pdevice=pdevice,
            )
        )
        names.append(
            self.add_plane(
                "blayer{}_zmin".format(ilayer),
                0,
                0,
                overt - dvert / 2,
                0,
                0,
                -1,
                pdevice=pdevice,
            )
        )
        names.append(
            self.add_plane(
                "blayer{}_zmax".format(ilayer),
                0,
                0,
                overt + dvert / 2,
                0,
                0,
                1,
                pdevice=pdevice,
            )
        )
        return ilayer, names


def materialname(material):
    try:
        name = material.name
    except AttributeError:
        name = material
    try:
        if "vacuum" in name.lower():
            name = "Vacuum"
    except AttributeError:
        name = "Vacuum"
    name = name.replace(" ", "_")
    return name


class Compositions(XrmcDevice):

    TYPE = "composition"

    def __init__(self, parent, name, atmosphere=None):
        super(Compositions, self).__init__(parent, name)
        self.clear()
        self.atmosphere = atmosphere

    @property
    def header(self):
        return [
            "Composition file",
            "",
            "List of compounds, defined by their density and elemental composition (weight fractions)",
        ]

    @property
    def code(self):
        if not self._dict:
            return []
        lines = super(Compositions, self).code
        for name, (elements, massfractions, density) in self._dict.items():
            lines += [
                None,
                ("Phase {}".format(name), "Name", True),
                ("NElem {}".format(len(elements)), " # elements", True),
            ]
            for formula, wfrac in zip(elements, massfractions):
                lines.append(
                    ("{} {}".format(formula, wfrac * 100), "mass fraction", True)
                )
            lines += [("Rho {}".format(density), "Mass density (g/cm3)", True)]
        return lines

    def add(self, name, elements, massfractions, density):
        self._dict[name] = list(map(str, elements)), massfractions, density

    def clear(self):
        self._dict = {}

    def add_material(self, material, name=None):
        """
        Args:
            material(str or Element/Compound/Mixture)
            name(str)
        Returns:
            str: xrmc composition name
        """
        if material is None:
            material = "vacuum"
        if name is None:
            name = materialname(material)
        if instance.isstring(material):
            if "vacuum" in material.lower():
                return "Vacuum"  # No need to add this
            material = compoundfromname(material)
        else:
            if "vacuum" in material.name.lower():
                return "Vacuum"  # No need to add this
        wfrac = material.elemental_massfractions()
        self.add(name, list(wfrac.keys()), list(wfrac.values()), material.density)
        return name

    @property
    def atmosphere(self):
        return self._atmosphere

    @atmosphere.setter
    def atmosphere(self, material):
        self._atmosphere = self.add_material(material)
        for objects in self.parent.objects(first=False):
            objects.change_atmosphere(material)


class XrmcWorldBuilder(object):
    def __init__(self, path, atmosphere=None):
        self.main = main = XrmcMain(path)
        self.compoundlib = main.add_device(
            Compositions, "compoundlib", atmosphere=atmosphere
        )
        self.quadrics = main.add_device(Quadrics, "quadric")
        self.objects = main.add_device(Objects, "objects")
        self.sample = main.add_device(Sample, "sample")

    def __repr__(self):
        return repr(self.main)

    def define_source(
        self, flux=None, energy=None, distance=None, beamsize=None, multiplicity=1
    ):
        self.main.remove_device(cls=Source)
        self.source = self.main.add_device(
            Source, "synchrotron", beamsize=beamsize, distance=distance
        )
        self.spectrum = self.main.add_device(
            Spectrum,
            "dcmspectrum",
            lines=[[energy, 0, flux]],
            multiplicity=multiplicity,
        )

    def add_diode(
        self,
        distance=None,
        activearea=None,
        ebinsize=None,
        polar=0,
        azimuth=0,
        poissonnoise=False,
        forcedetect=False,
        multiplicity=1,
        time=1,
    ):
        self.main.remove_device(cls=Detector)
        self.detector = self.main.add_device(
            SingleElementDetector,
            "detector",
            distance=distance,
            activearea=activearea,
            polar=polar,
            azimuth=azimuth,
            ebinsize=ebinsize,
            poissonnoise=poissonnoise,
            forcedetect=forcedetect,
            multiplicity=multiplicity,
            time=time,
        )

    def add_xrfdetector(
        self,
        distance=None,
        activearea=None,
        ebinsize=None,
        polar=0,
        azimuth=0,
        hoffset=0,
        voffset=0,
        poissonnoise=False,
        forcedetect=True,
        multiplicity=1,
        time=1,
        response=None,
        emin=None,
        emax=None,
        as_lines=None,
    ):
        if as_lines is None:
            as_lines = not bool(response)
        self.main.remove_device(cls=Detector)
        if response:
            cls = SDD
        else:
            cls = SingleElementDetector
            response = {}
        self.detector = self.main.add_device(
            cls,
            "detector",
            distance=distance,
            activearea=activearea,
            polar=polar,
            azimuth=azimuth,
            hoffset=hoffset,
            voffset=voffset,
            ebinsize=ebinsize,
            poissonnoise=poissonnoise,
            forcedetect=forcedetect,
            multiplicity=multiplicity,
            emin=emin,
            emax=emax,
            time=time,
            as_lines=as_lines,
            **response
        )

    def add_areadetector(
        self,
        distance=None,
        activearea=None,
        ebinsize=None,
        polar=0,
        azimuth=0,
        hpoffset=0,
        vpoffset=0,
        poissonnoise=False,
        forcedetect=True,
        multiplicity=1,
        time=1,
        pixelsize=None,
        dims=None,
    ):
        self.main.remove_device(cls=Detector)
        self.detector = self.main.add_device(
            AreaDetector,
            "detector",
            distance=distance,
            pixelsize=pixelsize,
            dims=dims,
            polar=polar,
            azimuth=azimuth,
            hpoffset=hpoffset,
            vpoffset=vpoffset,
            ebinsize=ebinsize,
            poissonnoise=poissonnoise,
            forcedetect=forcedetect,
            multiplicity=multiplicity,
            time=time,
        )

    def add_material(self, material, name=None):
        return self.main.compositions().add_material(material, name=name)

    def finalize(self, interactions=None):
        sample = self.sample
        sample.multiplicity = interactions
        self.main.save()

    def simulate(self, interactions=None):
        return self.main.simulate()

    # Legacy notebooks:
    def removesample(self):
        pass

    def addlayer(self, **kwargs):
        self.sample.add_layer(**kwargs)

    definesource = define_source
    adddiode = add_diode
    addxrfdetector = add_xrfdetector
    addareadetector = add_areadetector
    addmaterial = add_material

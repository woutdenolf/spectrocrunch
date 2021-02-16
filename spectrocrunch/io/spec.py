# -*- coding: utf-8 -*-

import re
import numpy as np
from collections import OrderedDict
import os
import logging

from ..utils import instance
from ..utils import units
from ..patch.pint import ureg
from .xiaedf import XiaNameParser
from ..process import axis

from PyMca5.PyMcaCore import SpecFileDataSource

logger = logging.getLogger(__name__)


def zapline_values(start, end, npixels):
    """Values of the pixel centers"""
    inc = (end - start) / np.float(npixels)
    return start + inc / 2 + inc * np.arange(npixels)


def zapline_range(start, end, npixels):
    """First and last pixel center"""
    inc = (end - start) / np.float(npixels)
    return [start + inc / 2, end - inc / 2]


def zapline_pixelsize(start, end, npixels):
    """Pixel size"""
    return (end - start) / np.float(npixels - 1)


def zapline_scansize(start, end, npixels):
    """Distance between last and first pixel center"""
    inc = (end - start) / np.float(npixels)
    return end - start - inc


def ascan_values(start, end, nsteps):
    """Values of the pixel centers"""
    inc = (end - start) / np.float(nsteps)
    return start + inc * np.arange(nsteps + 1)


def ascan_range(start, end, nsteps):
    """First and last pixel center"""
    return [start, end]


def ascan_pixelsize(start, end, nsteps):
    """Pixel size"""
    return (end - start) / np.float(nsteps)


def ascan_scansize(start, end, npixels):
    """Distance between last and first pixel center"""
    return end - start


def zapimage_submap(header, cmdlabel, scanrange, currentpos, microntounits):
    """
    Args:
        header(dict): edf header
        cmdlabel(str): scan command header key
        scanrange(num): in micron
        currentpos(dict): current motor positions
        microntounits(dict): factors to convert microns to motor units
    """

    # Old map motor values
    o = cmd_parser()
    result = o.parsezapimage(header[cmdlabel])
    if result["name"] != "zapimage":
        raise RuntimeError("Cannot extract zapimage information from edf header")
    fastvalues = zapline_values(
        result["startfast"], result["endfast"], result["npixelsfast"]
    )
    slowvalues = ascan_values(
        result["startslow"], result["endslow"], result["nstepsslow"]
    )

    # New map motor values
    pfasta = (
        currentpos[result["motfast"]]
        - scanrange * microntounits[result["motfast"]] / 2.0
    )
    pfastb = (
        currentpos[result["motfast"]]
        + scanrange * microntounits[result["motfast"]] / 2.0
    )
    pslowa = (
        currentpos[result["motslow"]]
        - scanrange * microntounits[result["motslow"]] / 2.0
    )
    pslowb = (
        currentpos[result["motslow"]]
        + scanrange * microntounits[result["motslow"]] / 2.0
    )

    ifasta = (np.abs(fastvalues - pfasta)).argmin()
    ifastb = (np.abs(fastvalues - pfastb)).argmin()
    islowa = (np.abs(slowvalues - pslowa)).argmin()
    islowb = (np.abs(slowvalues - pslowb)).argmin()

    result["startfast"] = fastvalues[ifasta]
    result["endfast"] = fastvalues[ifastb]
    result["npixelsfast"] = abs(ifastb - ifasta) + 1

    result["startslow"] = slowvalues[islowa]
    result["endslow"] = slowvalues[islowb]
    result["nstepsslow"] = abs(islowb - islowa)

    startpositions = {mot: header[mot] for mot in currentpos}
    startpositions[result["motfast"]] = result["startfast"]
    startpositions[result["motslow"]] = result["startslow"]

    d = (result["endfast"] - result["startfast"]) / (result["npixelsfast"] - 1.0)

    scancmd = "zapimage {} {} {} {} {} {} {} {} {} 0".format(
        result["motfast"],
        result["startfast"] + d / 2.0,
        result["endfast"] + d / 2.0,
        result["npixelsfast"],
        result["motslow"],
        result["startslow"],
        result["endslow"],
        result["nstepsslow"],
        int(result["time"].to("ms").magitude),
    )

    mvcmd = "mv " + " ".join(
        "{} {}".format(mot, pos) for mot, pos in startpositions.items()
    )

    if ifasta > ifastb:
        ifast = -1
    else:
        ifast = 1
    if islowa > islowb:
        islow = -1
    else:
        islow = 1

    for k, v in currentpos.items():
        if k != result["motfast"] and k != result["motslow"]:
            if v != startpositions[k]:
                logger.warning(
                    "Current position of {} ({}) is ignored and set to {}".format(
                        k, v, startpositions[k]
                    )
                )

    return scancmd, mvcmd, [[ifasta, ifastb + 1, ifast], [islowa, islowb + 1, islow]]


class cmd_parser(object):
    def __init__(self):
        self.fnumber = r"(?:[+-]?[0-9]*\.?[0-9]+)"
        self.inumber = r"\d+"
        self.blanks = r"\s+"
        self.motor = r"[a-zA-Z]+"
        self.motornum = r"[a-zA-Z0-9]+"

    def __call__(self, cmd):
        return self.parse(cmd)

    def parse(self, cmd):
        scanname = cmd.split(" ")[0]
        if scanname == "zapimage":
            return self.parsezapimage(cmd)
        elif scanname == "ascan":
            return self.parseascan(cmd)
        elif scanname == "zapline":
            return self.parsezapline(cmd)
        elif scanname == "zapenergy":
            return self.parsezapenergy(cmd)
        elif scanname == "mesh":
            return self.parsemesh(cmd)
        elif scanname == "puzzle":
            return self.parsepuzzle(cmd)
        else:
            return {"name": "unknown"}

    def match(self, cmd, patterns):
        for pattern in patterns:
            m = re.match(pattern, cmd)
            if m:
                return m
        return m

    def patternzapimage(self, name="zapimage"):
        pat1 = (
            "(?P<name>"
            + name
            + ")"
            + self.blanks
            + "(?P<motfast>"
            + self.motor
            + ")"
            + self.blanks
            + "(?P<startfast>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<endfast>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<npixelsfast>"
            + self.inumber
            + ")"
            + self.blanks
            + "(?P<time>"
            + self.inumber
            + ")"
            + self.blanks
            + "(?P<motslow>"
            + self.motor
            + ")"
            + self.blanks
            + "(?P<startslow>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<endslow>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<nstepsslow>"
            + self.inumber
            + ")"
        )

        pat2 = (
            "(?P<name>"
            + name
            + ")"
            + self.blanks
            + "(?P<motfast>"
            + self.motor
            + ")"
            + self.blanks
            + "(?P<startfast>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<endfast>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<npixelsfast>"
            + self.inumber
            + ")"
            + self.blanks
            + "(?P<motslow>"
            + self.motor
            + ")"
            + self.blanks
            + "(?P<startslow>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<endslow>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<nstepsslow>"
            + self.inumber
            + ")"
            + self.blanks
            + "(?P<time>"
            + self.inumber
            + ")"
        )

        return [pat1, pat2]

    def castzapimage(self, m):
        if m:
            result = m.groupdict()
            result["name"] = str(result["name"])

            result["motfast"] = str(result["motfast"])
            result["startfast"] = np.float(result["startfast"])
            result["endfast"] = np.float(result["endfast"])
            result["npixelsfast"] = np.int(result["npixelsfast"])

            result["motslow"] = str(result["motslow"])
            result["startslow"] = np.float(result["startslow"])
            result["endslow"] = np.float(result["endslow"])
            result["nstepsslow"] = np.int(result["nstepsslow"])

            result["time"] = ureg.Quantity(np.float(result["time"]), "ms")
        else:
            result = {"name": "unknown"}
        return result

    def patternpuzzle(self, name="puzzle"):
        return [
            "(?P<name>"
            + name
            + ")"
            + self.blanks
            + "(?P<motfast>"
            + self.motornum
            + ")"
            + self.blanks
            + "(?P<startfast>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<endfast>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<npixelsfast>"
            + self.inumber
            + ")"
            + self.blanks
            + "(?P<motslow>"
            + self.motornum
            + ")"
            + self.blanks
            + "(?P<startslow>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<endslow>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<nstepsslow>"
            + self.inumber
            + ")"
            + self.blanks
            + "(?P<time>"
            + self.inumber
            + ")"
        ]

    def castpuzzle(self, m):
        return self.castzapimage(m)

    def patternmesh(self, name="mesh"):
        return [
            "(?P<name>"
            + name
            + ")"
            + self.blanks
            + "(?P<motfast>"
            + self.motor
            + ")"
            + self.blanks
            + "(?P<startfast>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<endfast>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<nstepsfast>"
            + self.inumber
            + ")"
            + self.blanks
            + "(?P<motslow>"
            + self.motor
            + ")"
            + self.blanks
            + "(?P<startslow>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<endslow>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<nstepsslow>"
            + self.inumber
            + ")"
            + self.blanks
            + "(?P<time>"
            + self.fnumber
            + ")"
        ]

    def castmesh(self, m):
        if m:
            result = m.groupdict()
            result["name"] = str(result["name"])

            result["motfast"] = str(result["motfast"])
            result["startfast"] = np.float(result["startfast"])
            result["endfast"] = np.float(result["endfast"])
            result["nstepsfast"] = np.int(result["nstepsfast"])

            result["motslow"] = str(result["motslow"])
            result["startslow"] = np.float(result["startslow"])
            result["endslow"] = np.float(result["endslow"])
            result["nstepsslow"] = np.int(result["nstepsslow"])

            result["time"] = ureg.Quantity(np.float(result["time"]), "s")
        else:
            result = {"name": "unknown"}
        return result

    def patternzapenergy(self, name="zapenergy"):
        pat1 = (
            "(?P<name>"
            + name
            + ")"
            + self.blanks
            + "SUM"
            + self.blanks
            + "(?P<repeats>"
            + self.inumber
            + ")"
            + self.blanks
            + "(?P<time>"
            + self.fnumber
            + ")"
        )
        pat2 = (
            "(?P<name>"
            + name
            + ")"
            + self.blanks
            + "SUM2"
            + self.blanks
            + "(?P<repeats>"
            + self.inumber
            + ")"
            + self.blanks
            + "(?P<time>"
            + self.fnumber
            + ")"
        )
        return [pat1, pat2]

    def castzapenergy(self, m):
        if m:
            result = m.groupdict()
            result["name"] = str(result["name"])

            result["repeats"] = np.int(result["repeats"])
            result["time"] = ureg.Quantity(np.float(result["time"]), "ms")
        else:
            result = {"name": "unknown"}
        return result

    def patternzapline(self, name="zapline"):
        return [
            "(?P<name>"
            + name
            + ")"
            + self.blanks
            + "(?P<motfast>"
            + self.motor
            + ")"
            + self.blanks
            + "(?P<startfast>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<endfast>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<npixelsfast>"
            + self.inumber
            + ")"
            + self.blanks
            + "(?P<time>"
            + self.fnumber
            + ")"
        ]

    def castzapline(self, m):
        if m:
            result = m.groupdict()
            result["name"] = str(result["name"])

            result["motfast"] = str(result["motfast"])
            result["startfast"] = np.float(result["startfast"])
            result["endfast"] = np.float(result["endfast"])
            result["npixelsfast"] = np.int(result["npixelsfast"])

            result["time"] = ureg.Quantity(np.float(result["time"]), "ms")
        else:
            result = {"name": "unknown"}
        return result

    def patternascan(self, name="ascan"):
        return [
            "(?P<name>"
            + name
            + ")"
            + self.blanks
            + "(?P<motfast>"
            + self.motor
            + ")"
            + self.blanks
            + "(?P<startfast>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<endfast>"
            + self.fnumber
            + ")"
            + self.blanks
            + "(?P<nstepsfast>"
            + self.inumber
            + ")"
            + self.blanks
            + "(?P<time>"
            + self.fnumber
            + ")"
        ]

    def castascan(self, m):
        if m:
            result = m.groupdict()
            result["name"] = str(result["name"])

            result["motfast"] = str(result["motfast"])
            result["startfast"] = np.float(result["startfast"])
            result["endfast"] = np.float(result["endfast"])
            result["nstepsfast"] = np.int(result["nstepsfast"])

            result["time"] = ureg.Quantity(np.float(result["time"]), "s")
        else:
            result = {"name": "unknown"}
        return result

    def parsezapimage(self, cmd, name="zapimage"):
        patterns = self.patternzapimage(name=name)
        m = self.match(cmd, patterns)
        return self.castzapimage(m)

    def parsepuzzle(self, cmd, name="puzzle"):
        patterns = self.patternpuzzle(name=name)
        m = self.match(cmd, patterns)
        return self.castpuzzle(m)

    def parsemesh(self, cmd, name="mesh"):
        patterns = self.patternmesh(name=name)
        m = self.match(cmd, patterns)
        return self.castmesh(m)

    def parsezapenergy(self, cmd, name="zapenergy"):
        # Only SUM is called zapenergy, otherwise its called zapline
        patterns = self.patternzapenergy(name=name)
        m = self.match(cmd, patterns)
        return self.castzapenergy(m)

    def parsezapline(self, cmd, name="zapline"):
        patterns = self.patternzapline(name=name)
        m = self.match(cmd, patterns)
        return self.castzapline(m)

    def parseascan(self, cmd, name="ascan"):
        patterns = self.patternascan(name=name)
        m = self.match(cmd, patterns)
        return self.castascan(m)


class edfheader_parser(object):
    def __init__(
        self,
        fastlabel=None,
        slowlabel=None,
        speclabel=None,
        units=None,
        compensationmotors=None,
        axesnamemap=None,
        **otherlabels
    ):
        """
        Args:
            fastlabel(Optional(str)): "slow"
            slowlabel(Optional(str)): "fast"
            speclabel(Optional(str)): "Title"
            units(Optional(dict)): {'samz':'mm',...}
            compensationmotors(Optional(dict)): {'samz':['sampz'],...}
            axesnamemap(Optional(dict)): {'samz':'z',...}
        """
        self.fastlabel = fastlabel
        self.slowlabel = slowlabel
        self.speclabel = speclabel
        self.otherlabels = otherlabels
        if compensationmotors is None:
            compensationmotors = {}
        self.compensationmotors = compensationmotors
        if axesnamemap is None:
            axesnamemap = {}
        self.axesnamemap = axesnamemap
        self.specparser = cmd_parser()
        if units:
            self.units = units
        else:
            self.units = {}

    def __call__(self, cmd):
        return self.parse(cmd)

    def parse(self, header, defaultdims=None):
        out = {}
        if self.speclabel:
            try:
                out = self.specparser.parse(str(header[self.speclabel]))
            except KeyError:
                out = {}

        if self.fastlabel:
            try:
                out["motfast"] = str(header[self.fastlabel + "_mot"])
                out["startfast"] = np.float(header[self.fastlabel + "_start"])
                out["endfast"] = np.float(header[self.fastlabel + "_end"])
                out["npixelsfast"] = np.int(header[self.fastlabel + "_nbp"])
            except KeyError:
                pass

        if self.slowlabel:
            try:
                out["motslow"] = str(header[self.slowlabel + "_mot"])
                out["startslow"] = np.float(header[self.slowlabel + "_start"])
                out["endslow"] = np.float(header[self.slowlabel + "_end"])
                out["nstepsslow"] = np.int(header[self.slowlabel + "_nbp"])
            except KeyError:
                pass

        self._extract_floats(header, self.otherlabels, out)
        for k, labels in self.compensationmotors.items():
            for v in labels:
                self._extract_floats(header, {v: v}, out)

        # EDF: row-first ordering
        # data.shape == (Dim_2, Dim_1) == (slow, fast)
        axes = []
        if not defaultdims:
            defaultdims = (None, None)
        axes.append(
            self._extract_axis(
                header, out, "Dim_2", ndefault=defaultdims[0], fast=False
            )
        )
        axes.append(
            self._extract_axis(header, out, "Dim_1", ndefault=defaultdims[1], fast=True)
        )
        out["axes"] = axes

        if "name" not in out:
            out["name"] = "unknown"
            if len(axes) == 2:
                out["name"] = "zapimage"

        return out

    def _extract_floats(self, header, labels, out):
        for k, label in labels.items():
            if label:
                try:
                    u = self.units.get(label, "dimensionless")
                    out[k] = units.Quantity(np.float(header[label]), u)
                except KeyError:
                    if k not in out:
                        out[k] = ureg.Quantity(np.nan, u)

    def _extract_axis(self, header, out, dimkey, ndefault=None, fast=True):
        if fast:
            label = "fast"
            nlabel = "npixels"
        else:
            label = "slow"
            nlabel = "nsteps"
        try:
            name = out.get("mot" + label)
            u = self.units.get(name, "dimensionless")
            start = units.Quantity(out.pop("start" + label), units=u)
            end = units.Quantity(out.pop("end" + label), units=u)
            n = out.pop(nlabel + label)
            lst = self.compensationmotors.get(name, [])
            for mot in lst:
                pos = out.get(mot, None)
                if pos is not None:
                    if not np.isnan(pos.magnitude):
                        start += pos
                        end += pos
                        name = self.axesnamemap.get(name, name)
            if fast:
                return axis.zapscan(start, end, n, name=name)
            else:
                return axis.ascan(start, end, n, name=name)
        except KeyError:
            name = label
            if not ndefault:
                ndefault = int(header[dimkey])
            nsteps = ndefault - 1
            if fast:
                return axis.AxisRegular(0.5, nsteps + 0.5, nsteps, name=name)
            else:
                return axis.AxisRegular(0, nsteps, nsteps, name=name)


class spec(SpecFileDataSource.SpecFileDataSource):
    """An interface to a spec file"""

    def __init__(self, filename):
        """
        Args:
            filename(str): file name

        Raises:
            ValueError: file cannot be loaded
        """
        SpecFileDataSource.SpecFileDataSource.__init__(self, filename)
        self.parser = cmd_parser()

    def _parse_labels(self, labels):
        if instance.isstring(labels):
            return labels.lower()
        else:
            return [self._parse_labels(label) for label in labels]

    def _get_scan(self, scannumber):
        try:
            return self.getDataObject("{:d}.1".format(scannumber))
        except:
            msg = "Failed to retrieve scan number {} from {}".format(
                scannumber, self.sourceName
            )
            raise KeyError(msg)

    def _get_scan_info(self, scannumber):
        try:
            return self.getKeyInfo("{:d}.1".format(scannumber))
        except:
            msg = "Failed to retrieve scan number {} from {}".format(
                scannumber, self.sourceName
            )
            raise KeyError(msg)

    def _data_from_scan(self, scan, labels):
        lst = self._parse_labels(scan.info["LabelNames"])
        if not labels:
            return scan.data
        labels = self._parse_labels(labels)
        ind = [lst.index(label) for label in labels]
        if ind:
            return scan.data[:, ind]
        else:
            return None

    def labels(self, scannumber):
        return self._get_scan(scannumber).info["LabelNames"]

    @classmethod
    def _getxialoc(cls, header):
        loc = {
            "DIRECTORY": "",
            "RADIX": "",
            "ZAP SCAN NUMBER": "",
            "ZAP IMAGE NUMBER": "",
        }
        for s in header:
            if s.startswith("#C "):
                tmp = s[2:].split(":")
                if len(tmp) == 2:
                    tmp = [s.strip() for s in tmp]
                    if tmp[0] in loc:
                        loc[tmp[0]] = tmp[1]
            if s.startswith("#@XIAFILE "):
                filename = s[10:]
                loc["DIRECTORY"] = os.path.dirname(filename)
                filename = os.path.basename(filename).replace("X", "0")
                o = XiaNameParser().parse(filename)
                loc["RADIX"] = o.radix
                loc["ZAP SCAN NUMBER"] = o.mapnum
                loc["ZAP IMAGE NUMBER"] = 0
        return loc

    def getdata(self, scannumber, labelnames=None):
        """Get counters + info on saved data

        Args:
            scannumber(num): spec scan number
            labelnames(list(str)): list of counter names

        Returns:
            (np.array, dict): counters (nCounter x npts), info on collected data
        """
        scan = self._get_scan(scannumber)
        info = self._getxialoc(scan.info["Header"])
        data = self._data_from_scan(scan, labelnames)
        return data, info

    def getdata2(self, scannumber, labelnames=None):
        """Get counters

        Args:
            scannumber(num): spec scan number
            labelnames(list(str)): list of labels

        Returns:
            np.array: counters (nCounter x npts)
        """
        scan = self._get_scan(scannumber)
        return self._data_from_scan(scan, labelnames)

    def haslabels(self, scannumber, labels):
        scan = self._get_scan(scannumber)
        lst = self._parse_labels(scan.info["LabelNames"])
        labels = self._parse_labels(labels)
        if instance.isstring(labels):
            return labels in lst
        else:
            return [label in lst for label in labels]

    def getmotorvalues(self, scannumber, motors):
        """Get start positions for the specified motors"""
        info = self._get_scan_info(scannumber)
        motors = self._parse_labels(motors)
        names = self._parse_labels(info["MotorNames"])
        values = info["MotorValues"]
        if instance.isstring(motors):
            return values[names.index(motors)]
        else:
            return [values[names.index(mot)] for mot in motors]

    def hasmotors(self, scannumber, motors):
        info = self._get_scan_info(scannumber)
        motors = self._parse_labels(motors)
        names = self._parse_labels(info["MotorNames"])
        if instance.isstring(motors):
            return motors in names
        else:
            return [mot in names for mot in motors]

    def getxialocation(self, scannumber):
        """Get info on saved data"""
        info = self._get_scan_info(scannumber)
        ret = self._getxialoc(info["Header"])
        return ret

    def scancommand(self, scannumber):
        info = self._get_scan_info(scannumber)
        return self.parser.parse(info["Command"])

    def getdimensions(self, scannumber, motors):
        """Get scan dimensions for the specified motors"""
        info = self._get_scan_info(scannumber)
        motors = self._parse_labels(motors)
        names = self._parse_labels(info["MotorNames"])
        values = info["MotorValues"]
        cmd = info["Command"]

        # Parse motor positions
        ret = {mot: values[names.index(mot)] for mot in motors}

        # Parse command
        result = self.parser.parse(cmd)
        if result["name"] == "zapimage":
            if result["motfast"] in motors:
                ret[result["motfast"]] = np.array(
                    zapline_range(
                        result["startfast"], result["endfast"], result["npixelsfast"]
                    )
                )
            if result["motslow"] in motors:
                ret[result["motslow"]] = np.array(
                    ascan_range(
                        result["startslow"], result["endslow"], result["nstepsslow"]
                    )
                )
        elif result["name"] == "ascan":
            if result["motfast"] in motors:
                ret[result["motfast"]] = np.array(
                    ascan_range(
                        result["startfast"], result["endsfast"], result["nstepsfast"]
                    )
                )
        elif result["name"] == "mesh":
            if result["motfast"] in motors:
                ret[result["motfast"]] = np.array(
                    ascan_range(
                        result["startfast"], result["endfast"], result["npixelsfast"]
                    )
                )
            if result["motslow"] in motors:
                ret[result["motslow"]] = np.array(
                    ascan_range(
                        result["startslow"], result["endslow"], result["nstepsslow"]
                    )
                )

        return ret

    @staticmethod
    def addunit(value, mot, units):
        u = units.get(mot, None)
        if units is None:
            return value
        else:
            return ureg.Quantity(value, u)

    def getimages(self, motors=None, units={}):
        """Get list of all images"""
        ret = OrderedDict()
        for k in self.getSourceInfo()["KeyList"]:

            info = self.getKeyInfo(k)
            cmd = info["Command"]
            # or cmd.startswith("mesh"): NO LINK TO DATA FOR MESH
            if cmd.startswith("zapimage") or cmd.startswith("puzzle"):
                result = self.parser.parse(cmd)
                if result["name"] == "unknown":
                    continue

                add = {}
                add["specnumber"] = k.split(".")[0]
                add["scansize"] = "{} x {}".format(
                    result["npixelsfast"], result["nstepsslow"] + 1
                )

                def ffast(op):
                    return self.addunit(
                        op(
                            result["startfast"],
                            result["endfast"],
                            result["npixelsfast"],
                        ),
                        result["motfast"],
                        units,
                    )

                def fslow(op):
                    return self.addunit(
                        op(
                            result["startslow"], result["endslow"], result["nstepsslow"]
                        ),
                        result["motslow"],
                        units,
                    )

                ufast = units.get(result["motfast"], None)
                uslow = units.get(result["motslow"], None)

                add["scansize_units"] = "{} x {}".format(
                    ffast(zapline_scansize), ffast(ascan_scansize)
                )

                add["pixelsize"] = "{} x {}".format(
                    ffast(zapline_pixelsize), ffast(ascan_pixelsize)
                )

                add["puzzle"] = result["name"] == "puzzle"
                add["mesh"] = result["name"] == "mesh"

                if motors is None:
                    add["motors"] = []
                else:
                    motors = self._parse_labels(motors)
                    names = self._parse_labels(info["MotorNames"])
                    values = info["MotorValues"]
                    add["motors"] = [
                        values[names.index(mot)] if mot in names else np.nan
                        for mot in motors
                    ]

                if result["name"] != "mesh":
                    h = info["Header"]
                    h = [i.split(":")[1].strip() for i in h if "#C " in i]
                    add["scandir"] = h[0]
                    add["scanname"] = h[1]
                    add["scannumber"] = h[2]

                    ret["{}_{}".format(h[1], int(h[2]))] = add

        return ret

    def getregexscans(self, pattern):
        ret = OrderedDict()
        fmt = re.compile(pattern)
        for k in self.getSourceInfo()["KeyList"]:
            info = self.getKeyInfo(k)
            m = fmt.match(info["Command"])
            if m:
                ret[int(k.split(".")[0])] = m.groupdict()
        return ret

    def extractxanesinfo(self, skip=None, nrmin=None, nrmax=None):
        """Get list of all ID21 XANES"""
        ret = []

        lst = self.getSourceInfo()["KeyList"]
        if skip is None:
            skip = []

        for k in lst:
            scannumber = int(k.split(".")[0])
            if nrmin is not None:
                if scannumber < nrmin:
                    continue
            if nrmax is not None:
                if scannumber > nrmax:
                    continue
            if scannumber in skip:
                continue

            info = self.getKeyInfo(k)
            if info["Command"].startswith("zapline mono "):
                if info["Lines"] == 0:
                    continue
                ret += [{"scannumber": scannumber, "repeats": 1}]

            elif info["Command"].startswith("zapenergy SUM "):
                if info["Lines"] == 0:
                    continue
                ret += [
                    {
                        "scannumber": scannumber,
                        "repeats": int(info["Command"].split(" ")[2]),
                    }
                ]
            elif info["Command"].startswith("zapenergy SUM2 "):
                # This used to be the sum fo the repeats, energy interpolated
                if info["Lines"] == 0:
                    continue
                ret[-1] = {
                    "scannumber": scannumber,
                    "repeats": int(info["Command"].split(" ")[2]),
                }

        return ret

    def extractxanes(self, scannumbers, labelnames):
        """Get list of specific ID21 XANES with counters"""
        ret = {}

        for scannumber in scannumbers:
            k = "{:d}.1".format(scannumber)

            info = self.getKeyInfo(k)
            if info["Command"].startswith("zapline mono"):
                result = self.parser.parse(info["Command"])
                if result["name"] != "zapline" or info["Lines"] == 0:
                    continue
                ret[scannumber] = {
                    "repeats": 1,
                    "data": self.getdata2(scannumber, labelnames),
                    "time": result["time"],
                    "labels": ["{}.{}".format(s, scannumber) for s in labelnames],
                }
            elif info["Command"].startswith("zapenergy SUM"):
                result = self.parser.parse(info["Command"])
                if result["name"] != "zapenergy" or info["Lines"] == 0:
                    continue
                ret[scannumber] = {
                    "repeats": int(info["Command"].split(" ")[2]),
                    "data": self.getdata2(scannumber, labelnames),
                    "time": result["time"],
                    "labels": ["{}.{}".format(s, scannumber) for s in labelnames],
                }

        return ret

    def extractxanesginfo(
        self,
        keepsum=False,
        sumingroups=False,
        keepindividual=False,
        skip=None,
        nrmin=None,
        nrmax=None,
    ):
        """Get list of all ID21 XANES, grouping repeats"""
        data = self.extractxanesinfo(skip=skip, nrmin=nrmin, nrmax=nrmax)

        ret = []

        bproc = [True] * len(data)

        # Groups: [rep1,rep2,...,(sum)]
        for i in range(len(data)):
            if not bproc[i]:
                continue

            n = data[i]["repeats"]
            if n > 1:
                # [rep1,rep2,....]
                i0 = i
                i1 = i
                while (data[i0 - 1]["repeats"] if i0 > 0 else 0) == 1:
                    i0 -= 1
                rng = range(max(i0, i - n), i1)
                add = [data[k]["scannumber"] for k in rng if bproc[k]]
                for l in rng:
                    bproc[l] = keepindividual

                # [rep1,rep2,...,(sum)]
                if keepsum:
                    if sumingroups:
                        add += [data[i1]["scannumber"]]
                        bproc[i1] = keepindividual
                else:
                    bproc[i1] = False

                if add:
                    ret += [add]

        # Add each scan as an individual group
        for i in range(len(data)):
            if bproc[i]:
                ret += [[data[i]["scannumber"]]]

        return ret

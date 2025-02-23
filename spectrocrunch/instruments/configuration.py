import os
from ..patch.pint import ureg
from ..utils.classfactory import with_metaclass
from ..utils import listtools


# defaultdict causes jsonpickle issues
class UnitDict(dict):
    def __getitem__(self, index):
        try:
            return super(UnitDict, self).__getitem__(index)
        except KeyError:
            return ureg.dimensionless


class InstrumentInfo(with_metaclass(object)):
    def __init__(self, **info):
        data = {}
        data["imagemotors"] = info.get("imagemotors", [])
        data["imageaxes"] = info.get("imageaxes", ("y", "x"))
        units = UnitDict()
        units.update(info.get("units", {}))
        data["units"] = units
        data["encoderinfo"] = info.get("encoderinfo", {})  # steps/motor units
        data["compensationmotors"] = info.get("compensationmotors", {})
        edfheaderkeys = {
            "speclabel": "title",
            "energylabel": "energy",
            "timelabel": "time",
            "energyunit": "keV",
            "timeunit": "s",
            "fastlabel": "fast",
            "slowlabel": "slow",
        }
        edfheaderkeys.update(info.get("edfheaderkeys", {}))
        data["edfheaderkeys"] = edfheaderkeys
        counterdict = {"I0_counts": "I0", "It_counts": "It"}
        counterdict.update(info.get("counterdict", []))
        data["counterdict"] = counterdict
        data["counter_reldir"] = info.get("counter_reldir", ".")
        data["metadata"] = info.get("metadata", "counters")
        data["xraysource"] = info.get("xraysource", "synchrotron")
        data["_diodeI0"] = info.get("diodeI0", {})
        data["_diodeIt"] = info.get("diodeIt", {})
        data["_optics"] = info.get("optics", {})
        data["_specmotornames"] = info.get("specmotornames", {})
        data["speccounternames"] = info.get("speccounternames", {})
        data["h5stackgroups"] = info.get(
            "h5stackgroups", ["counters", r"^detector(\d+|sum)$"]
        )
        datalocationinfo = info.get("datalocation", {})
        for k in ["xrf_dynamic", "xrf_static", "ff_dynamic", "ff_static"]:
            if k not in datalocationinfo:
                datalocationinfo[k] = {"subdir": ""}
        data["datalocationinfo"] = datalocationinfo
        self._data = data

    def __getattr__(self, key):
        if key in self._data:
            return self._data[key]
        else:
            raise AttributeError(key)

    def __getstate__(self):
        return {"_data": self._data}

    def __setstate__(self, state):
        self._data = state["_data"]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self._data == other._data
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def counters(self, include=None, exclude=None):
        if include is None:
            lst = self.counterdict.keys()
        else:
            lst = include
        if exclude is None:
            exclude = []
        lst = [k for k in lst if k not in exclude]
        ret = [self.counterdict[k] for k in lst if k in self.counterdict]
        return list(listtools.flatten(ret))

    @property
    def metadata(self):
        return self._metadata

    @metadata.setter
    def metadata(self, value):
        if value == "xia":
            self._metadata = "xia"
        else:
            self._metadata = "counters"

    def _devicename(self, device, position_dict):
        """
        Args:
            device(str or None): empty string: no device used
                                 None: don't know if/which device is used
            position_dict(dict): motor positions
        Returns:
            str
        """
        name = None
        for motname, limits in device.items():
            pos = position_dict.get(motname, None)
            if pos is not None:
                namea, posa, nameb, posb = limits
                if abs(pos - posa) < abs(pos - posb):
                    name = namea
                else:
                    name = nameb
                if name:
                    return name
        if "default" in device:
            name = device["default"]
        return name

    def diodeI0(self, position_dict):
        return self._devicename(self._diodeI0, position_dict)

    def diodeIt(self, position_dict):
        return self._devicename(self._diodeIt, position_dict)

    def optics(self, position_dict):
        return self._devicename(self._optics, position_dict)

    def specmotornames(self, names):
        return [self._specmotornames[k] for k in names]

    def motornames(self):
        return self._specmotornames.keys()

    def xrfsubdir(self, type="dynamic"):
        return self.datalocationinfo["xrf_{}".format(type)]["subdir"]

    def ffsubdir(self, type="dynamic"):
        return self.datalocationinfo["ff_{}".format(type)]["subdir"]

    def xrflocation(self, sample, dataset, **kwargs):
        return self.datalocation(sample, dataset, self.xrfsubdir(**kwargs))

    def fflocation(self, sample, dataset, **kwargs):
        return self.datalocation(sample, dataset, self.ffsubdir(**kwargs))

    def datalocation(self, sample, dataset, subdir):
        # root/sample/sample_dataset/subdir/sample_dataset*.*
        radix = "_".join([k for k in [sample, dataset] if k])
        subdir = os.path.join(sample, radix, subdir)
        return radix, subdir

    def edfparser_info(self):
        return {
            "time": self.edfheaderkeys.get("timelabel", None),
            "energy": self.edfheaderkeys.get("energylabel", None),
            "speclabel": self.edfheaderkeys.get("speclabel", None),
            "slowlabel": self.edfheaderkeys.get("slowlabel", None),
            "fastlabel": self.edfheaderkeys.get("fastlabel", None),
            "units": self.units,
            "compensationmotors": self.compensationmotors,
            "axesnamemap": self.imagemotors,
        }


class ESRF_ID21_SXM(InstrumentInfo):
    longname = "ESRF ID21: Scanning X-ray Microscope"
    shortname = "id21"
    aliases = ["sxm", "id21", "ID21"]

    def __init__(self, **info):
        info["imagemotors"] = info.get(
            "imagemotors", {"samy": "y", "sampy": "y", "samz": "z", "sampz": "z"}
        )
        info["imageaxes"] = ("z", "y")
        info["compensationmotors"] = info.get(
            "compensationmotors",
            {
                "samy": ["sampy"],
                "samz": ["sampz"],
                "sampy": ["samy"],
                "sampz": ["samz"],
            },
        )
        info["encoderinfo"] = info.get(
            "encoderinfo",
            {
                "samy": {"resolution": 52500, "counter": "arr_samy"},
                "samz": {"resolution": 50000, "counter": "arr_samz"},
            },
        )

        info["edfheaderkeys"] = info.get(
            "edfheaderkeys", {"speclabel": "Title", "energylabel": "DCM_Energy"}
        )

        info["counterdict"] = info.get(
            "counterdict",
            {
                "I0_counts": "arr_iodet",
                "It_counts": "arr_idet",
                "If_counts": "arr_fdet",
                "xrficr": "xmap_icr",
                "xrfocr": "xmap_ocr",
                "xrfdt": "xmap_dt",
                "encoders": ["arr_samy", "arr_samz"],
                "xrfroi": [
                    "xmap_x1",
                    "xmap_x1c",
                    "xmap_x2",
                    "xmap_x2c",
                    "xmap_x3",
                    "xmap_x3c",
                ],
                "calc": ["arr_absorp1", "arr_absorp2", "arr_absorp3"],
                "counters": ["arr_", "xmap_"],
            },
        )

        info["diodeI0"] = info.get(
            "diodeI0",
            {"istopz": ("iodet2", -20, "", -1.3), "ioz": ("iodet1", 7, "", 23)},
        )

        info["diodeIt"] = info.get("diodeIt", {"detz": ("idet", 10, "", 21)})

        info["optics"] = info.get("optics", {"zpz": ("KB", 7, "", 6.5)})

        info["specmotornames"] = info.get(
            "specmotornames",
            {
                "ioz": "diodeIoZ",
                "istopz": "istopz",
                "energy": "Energy MONO",
                "zpz": "zpz",
                "detz": "detz",
            },
        )

        info["speccounternames"] = info.get(
            "speccounternames",
            {
                "I0_counts": "iodet",
                "It_counts": "idet",
                "energy": "energym",
                "time": "seconds",
                "It_photons": "photons",
            },
        )

        info["units"] = info.get("units", {})
        info["units"]["time"] = info["units"].get("time", ureg.seconds)
        info["units"]["energy"] = info["units"].get("energy", ureg.Unit("keV"))
        info["units"]["DCM_Energy"] = info["units"].get("DCM_Energy", ureg.Unit("keV"))
        info["units"]["I0_counts"] = info["units"].get("I0_counts", ureg.dimensionless)
        info["units"]["It_counts"] = info["units"].get("It_counts", ureg.dimensionless)
        info["units"]["I0_cps"] = info["units"].get("I0_cps", ureg.Unit("Hz"))
        info["units"]["It_cps"] = info["units"].get("It_cps", ureg.Unit("Hz"))
        info["units"]["I0_current"] = info["units"].get("I0_current", ureg.Unit("Hz"))
        info["units"]["It_current"] = info["units"].get("It_current", ureg.Unit("Hz"))
        info["units"]["I0_photons"] = info["units"].get(
            "I0_photons", ureg.dimensionless
        )
        info["units"]["It_photons"] = info["units"].get(
            "It_photons", ureg.dimensionless
        )
        info["units"]["It_flux"] = info["units"].get("It_flux", ureg.Unit("Hz"))
        info["units"]["I0_flux"] = info["units"].get("I0_flux", ureg.Unit("Hz"))
        info["units"]["samy"] = info["units"].get("samy", ureg.millimeter)
        info["units"]["samz"] = info["units"].get("samz", ureg.millimeter)
        info["units"]["samx"] = info["units"].get("samx", ureg.millimeter)
        info["units"]["sampz"] = info["units"].get("sampz", ureg.micrometer)
        info["units"]["sampy"] = info["units"].get("sampy", ureg.micrometer)

        info["datalocation"] = info.get(
            "datalocation",
            {
                "xrf_dynamic": {"subdir": "zap"},
                "xrf_static": {"subdir": "xrf"},
                "ff_dynamic": {"subdir": "ff"},
                "ff_static": {"subdir": "ff"},
            },
        )

        super(ESRF_ID21_SXM, self).__init__(**info)


class ESRF_ID21_MICRODIFF(InstrumentInfo):
    longname = "ESRF ID21: Micro-Diffraction end-station"
    shortname = "id21"
    aliases = ["microdiff"]

    def __init__(self, **info):
        info["imagemotors"] = info.get(
            "imagemotors", {"samh": "h", "samph": "h", "samv": "v", "sampv": "v"}
        )
        info["imageaxes"] = ("v", "h")
        info["units"] = info.get(
            "units",
            {
                "samh": ureg.millimeter,
                "samv": ureg.millimeter,
                "samd": ureg.millimeter,
                "samph": ureg.micrometer,
                "sampv": ureg.micrometer,
            },
        )
        info["compensationmotors"] = info.get(
            "compensationmotors",
            {
                "samh": ["samph"],
                "samv": ["sampv"],
                "samph": ["samh"],
                "sampv": ["samv"],
            },
        )

        info["counterdict"] = info.get(
            "counterdict",
            {
                "I0_counts": "zap_iodet",
                "It_counts": "zap_idet",
                "xrficr": "xmap_icr",
                "xrfocr": "xmap_ocr",
                "xrfdt": "xmap_dt",
                "xrfroi": [
                    "xmap_x1",
                    "xmap_x1c",
                    "xmap_x2",
                    "xmap_x2c",
                    "xmap_x3",
                    "xmap_x3c",
                ],
            },
        )

        info["datalocation"] = info.get(
            "datalocation",
            {
                "xrf_dynamic": {"subdir": "zap"},
                "xrf_static": {"subdir": "xrf"},
                "xrd_dynamic": {"subdir": "zap"},
                "xrd_static": {"subdir": "xrd"},
            },
        )

        super(ESRF_ID21_MICRODIFF, self).__init__(**info)


class ESRF_ID16B(InstrumentInfo):
    longname = "ESRF ID16B"
    shortname = "id16b"
    aliases = ["id16b", "ID16b", "ID16B", "ID16NA"]

    def __init__(self, **info):
        info["imagemotors"] = info.get(
            "imagemotors", {"sy": "y", "sz": "z", "sampy": "y", "sampz": "z"}
        )
        info["imageaxes"] = ("z", "y")
        info["compensationmotors"] = info.get(
            "compensationmotors",
            {"sy": ["sampy"], "sz": ["sampz"], "sampy": ["sy"], "sampz": ["sz"]},
        )

        info["diodeI0"] = info.get(
            "diodeI0", {"default": "id16b_IC", "unknown": "id16b_I0"}
        )

        info["edfheaderkeys"] = info.get(
            "edfheaderkeys", {"speclabel": "title", "energylabel": "energy"}
        )

        info["diodeIt"] = info.get("diodeIt", {"default": "id16b_It"})

        info["optics"] = info.get("optics", {"zpz": "KB"})

        info["counterdict"] = info.get(
            "counterdict",
            {
                "I0_counts": "zap_p201_I0",
                "It_counts": "zap_p201_It",
                "IC_counts": "zap_p201_IC",
                "counters": ["zap_", "xmap_"],
            },
        )

        info["speccounternames"] = info.get(
            "speccounternames",
            {
                "IC_counts": "p201_IC",
                "I0_counts": "p201_I0",
                "It_counts": "p201_It",
                "IC_current": "IC",
                "I0_current": "I0",
                "It_current": "It",
                "energy": "energy",
                "time": "Seconds",
                "It_flux": "flux_It",
            },
        )

        info["units"] = info.get("units", {})
        info["units"]["time"] = info["units"].get("time", ureg.seconds)
        info["units"]["I0_counts"] = info["units"].get("I0_counts", ureg.dimensionless)
        info["units"]["It_counts"] = info["units"].get("It_counts", ureg.dimensionless)
        info["units"]["IC_counts"] = info["units"].get("IC_counts", ureg.dimensionless)
        info["units"]["I0_cps"] = info["units"].get("I0_cps", ureg.Unit("Hz"))
        info["units"]["It_cps"] = info["units"].get("It_cps", ureg.Unit("Hz"))
        info["units"]["IC_cps"] = info["units"].get("IC_cps", ureg.Unit("Hz"))
        info["units"]["energy"] = info["units"].get("energy", ureg.Unit("keV"))
        info["units"]["It_flux"] = info["units"].get("It_flux", ureg.Unit("Hz"))
        info["units"]["I0_flux"] = info["units"].get("I0_flux", ureg.Unit("Hz"))
        info["units"]["I0_current"] = info["units"].get("I0_current", ureg.Unit("pA"))
        info["units"]["It_current"] = info["units"].get("It_current", ureg.Unit("pA"))
        info["units"]["IC_current"] = info["units"].get("IC_current", ureg.Unit("pA"))
        info["units"]["sx"] = info["units"].get("sx", ureg.millimeter)
        info["units"]["sy"] = info["units"].get("sy", ureg.millimeter)
        info["units"]["sz"] = info["units"].get("sz", ureg.millimeter)
        info["units"]["sampz"] = info["units"].get("sampz", ureg.micrometer)
        info["units"]["sampy"] = info["units"].get("sampy", ureg.micrometer)

        info["metadata"] = info.get("metadata", "xia")

        super(ESRF_ID16B, self).__init__(**info)

    def datalocation(self, sample, dataset, subdir):
        # root/sample/subdir/dataset*.*
        if not dataset:
            dataset = sample
        radix = dataset
        subdir = os.path.join(sample, subdir)
        return radix, subdir


factory = InstrumentInfo.factory
registry = InstrumentInfo.clsregistry


def getinstrument(instrument=None, instrument_parameters=None, **kwargs):
    if instrument is None:
        raise RuntimeError("You need to specify an instrument.")
    if isinstance(instrument, InstrumentInfo):
        return instrument
    else:
        if instrument_parameters is None:
            instrument_parameters = {}
        return factory(instrument, **instrument_parameters)

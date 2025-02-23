import os
import logging

from . import xrf as xrfgeometries
from ..detectors import diode as diodes
from ..detectors import xrf as xrfdetectors
from ..sources import xray as xraysources
from ..optics import xray as xrayoptics
from ..instruments import configuration
from ..utils import units
from ..utils import instance
from ..utils.classfactory import with_metaclass
from ..io import spec
from ..math.linop import LinearOperator
from ..math.utils import weightedsum
from ..utils.copyable import Copyable

import numpy as np
import matplotlib.pyplot as plt
from pint import errors as pinterrors

logger = logging.getLogger(__name__)


def defaultunit():
    return units.dimensionless


class QXRFGeometry(with_metaclass(Copyable)):
    """Quantitative XRF geometry with I0 and It diodes"""

    def __init__(
        self,
        instrument=None,
        diodeI0=None,
        diodeIt=None,
        simplecalibration=True,
        fluxid=None,
        transmissionid=None,
        xrfgeometries=None,
        xrfdetector=None,
        xrfgeometry=None,
        optics=None,
        source=None,
        reference=units.Quantity(1e10, "hertz"),
        defaultreferencetime=None,
        defaultexpotime=units.Quantity(0.1, "s"),
    ):
        """
        Args:
            instrument(str)
            diodeI0(str)
            diodeIt(str)
            simplecalibration(bool)
            fluxid(str): prefix in diode calibration parameters
                         fluxid_counts, fluxid_flux, ...
            transmissionid(str): prefix in diode calibration parameters
                                 transmissionid_counts, transmissionid_flux, ...
            xrfgeometries(Optional(list))
            xrfdetector(Optional(str))
            xrfgeometry(Optional(str))
            optics(Optional(str))
            source(Optional(str)): synchrotron by default
            reference(Optional(Quantity)): normalize XRF data to this I0
            defaultreferencetime(Optional(Quantity)): normalize XRF data to this expo time
            defaultexpotime(Optional(Quantity)): XRF expo time in case not provided
        """
        self._prepare_update_devices()
        self.instrument = instrument
        self.simplecalibration = simplecalibration
        self.optics = optics
        if source is None:
            source = "synchrotron"
        self.source = source
        if xrfgeometries is None:
            if xrfgeometry is not None and xrfdetector is not None:
                xrfgeometries = [(xrfgeometry, xrfdetector)]
            else:
                xrfgeometries = []
        self.xrfgeometries = xrfgeometries
        self.diodeI0 = diodeI0
        self.diodeIt = diodeIt
        self.reference = reference  # could be flux or counts
        self.defaultreferencetime = defaultreferencetime
        self.defaultexpotime = defaultexpotime
        if fluxid is None:
            fluxid = "I0"
        self.fluxid = fluxid
        if transmissionid is None:
            transmissionid = "It"
        self.transmissionid = transmissionid

    def __getstate__(self):
        return {
            "instrument": self.instrument,
            "simplecalibration": self.simplecalibration,
            "optics": self.optics,
            "source": self.source,
            "xrfgeometries": self.xrfgeometries,
            "diodeI0": self.diodeI0,
            "diodeIt": self.diodeIt,
            "reference": self.reference,
            "defaultreferencetime": self.defaultreferencetime,
            "defaultexpotime": self.defaultexpotime,
            "fluxid": self.fluxid,
            "transmissionid": self.transmissionid,
        }

    def __setstate__(self, state):
        self._prepare_update_devices()
        self.instrument = state["instrument"]
        self.simplecalibration = state["simplecalibration"]
        self.optics = state["optics"]
        self.source = state["source"]
        self.xrfgeometries = state["xrfgeometries"]
        self.diodeI0 = state["diodeI0"]
        self.diodeIt = state["diodeIt"]
        self.reference = state["reference"]
        self.defaultreferencetime = state["defaultreferencetime"]
        self.defaultexpotime = state["defaultexpotime"]
        self.fluxid = state["fluxid"]
        self.transmissionid = state["transmissionid"]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (
                self.instrument == other.instrument
                and self.simplecalibration == other.simplecalibration
                and self.optics == other.optics
                and self.source == other.source
                and self.xrfgeometries == other.xrfgeometries
                and self.diodeI0 == other.diodeI0
                and self.diodeIt == other.diodeIt
                and self.reference == other.reference
                and self.defaultreferencetime == other.defaultreferencetime
                and self.defaultexpotime == other.defaultexpotime
                and self.fluxid == other.fluxid
                and self.transmissionid == other.transmissionid
            )
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        lst = []

        if self.diodeI0 is None:
            diodeI0 = None
        else:
            diodeI0 = self.diodeI0
            if hasattr(diodeI0, "geometry"):
                diodeI0 = diodeI0.geometry
        lst.append("DiodeI0:\n========\n{}".format(diodeI0))

        if self.diodeIt is None:
            diodeIt = None
        else:
            diodeIt = self.diodeIt
            if hasattr(diodeIt, "geometry"):
                diodeIt = diodeIt.geometry
        lst.append("DiodeIt:\n========\n{}".format(diodeIt))

        for i, geom in enumerate(self.xrfgeometries, 1):
            lst.append("XRF geometry ({}):\n=================\n{}".format(i, geom))

        lst.append("Flux monitor:\n===========")
        try:
            self.reference.to("Hz")
            refname = "Normalize to flux"
        except pinterrors.DimensionalityError:
            refname = "Normalize to diode counts"
        lst.append("{}:\n {:~e}".format(refname, self.reference))
        referencetime = self.defaultreferencetime
        if referencetime is None:
            referencetime = "<data>"
        else:
            referencetime = "{:~}".format(referencetime)
        lst.append("Normalize to exposure time:\n {}".format(referencetime))
        lst.append("Default exposure time:\n {:~}".format(self.defaultexpotime))

        return "\n".join(lst)

    @property
    def instrument(self):
        return self._instrument

    @instrument.setter
    def instrument(self, value):
        self._instrument = configuration.getinstrument(instrument=value)

    @property
    def simplecalibration(self):
        return self._simplecalibration

    @simplecalibration.setter
    def simplecalibration(self, value):
        self._simplecalibration = value
        self._update_devices()

    @property
    def diodeI0(self):
        return self._diodeI0

    @diodeI0.setter
    def diodeI0(self, device):
        self._diodeI0 = self._generate_device(device, diodes.factory)
        self._update_devices()

    @property
    def diodeIt(self):
        return self._diodeIt

    @diodeIt.setter
    def diodeIt(self, device):
        self._diodeIt = self._generate_device(device, diodes.factory)
        self._update_devices()

    @property
    def source(self):
        return self._source

    @source.setter
    def source(self, device):
        self._source = self._generate_device(device, xraysources.factory)
        self._update_devices()

    @property
    def optics(self):
        return self._optics

    @optics.setter
    def optics(self, device):
        self._optics = self._generate_device(device, xrayoptics.factory)
        self._update_devices()

    @property
    def xrfgeometries(self):
        return self._xrfgeometries

    @property
    def xrfgeometry(self):
        if len(self._xrfgeometries) != 1:
            raise RuntimeError("More than one XRF geometry, use 'xrfgeometries'")
        return self.xrfgeometries[0]

    @xrfgeometries.setter
    def xrfgeometries(self, values):
        """
        Args:
            values(list): list of geometry objects and/or
                          (geometry, detector) tuples.
                          The tuple elements can be of type str or dict.
                          In case of dict: {'name':'...', 'parameters':{...}}
        """
        self._xrfgeometries = []
        for geom in values:
            if not isinstance(geom, xrfgeometries.XRFGeometry):
                xrfgeometry, xrfdetector = geom

                # Instantiate geometry:
                if instance.isstring(xrfgeometry):
                    name = xrfgeometry
                    kwargs = {}
                else:
                    name = xrfgeometry["name"]
                    kwargs = xrfgeometry.get("parameters", {})
                geom = self._generate_device(name, xrfgeometries.factory, **kwargs)
                if not geom:
                    continue

                # Instantiate detector:
                if instance.isstring(xrfdetector):
                    name = xrfdetector
                    kwargs = {}
                else:
                    name = xrfdetector["name"]
                    kwargs = xrfdetector.get("parameters", {})
                geom.detector = self._generate_device(
                    name, xrfdetectors.factory, **kwargs
                )
            self._xrfgeometries.append(geom)
        self._update_devices()

    def _generate_device(self, device, factory, *args, **kwargs):
        if instance.isstring(device):
            return factory(device, *args, **kwargs)
        else:
            return device

    def _prepare_update_devices(self):
        self._diodeI0 = None
        self._diodeIt = None
        self._xrfgeometries = []

    def _update_devices(self):
        if self.diodeI0 is not None:
            self.diodeI0.optics = self.optics
            if self.diodeI0.secondarytarget is not None:
                self.diodeI0.secondarytarget.geometry.source = self.source
            self.diodeI0.simplecalibration = self.simplecalibration
        if self.xrfgeometries and self.source is not None:
            for geom in self.xrfgeometries:
                geom.source = self.source

    @property
    def xrf_positions(self):
        return units.asqarray([geom.detectorposition for geom in self.xrfgeometries])

    @xrf_positions.setter
    def xrf_positions(self, value):
        for geom, pos in zip(self.xrfgeometries, value):
            geom.detectorposition = pos

    @property
    def xrf_distances(self):
        return units.asqarray([geom.distance for geom in self.xrfgeometries])

    @xrf_distances.setter
    def xrf_distances(self, value):
        for geom, pos in zip(self.xrfgeometries, value):
            geom.distance = pos

    @property
    def xrf_activeareas(self):
        return units.asqarray([geom.detector.activearea for geom in self.xrfgeometries])

    @xrf_activeareas.setter
    def xrf_activeareas(self, value):
        for geom, pos in zip(self.xrfgeometries, value):
            geom.detector.activearea = pos

    @property
    def xrf_anglesin(self):
        return [geom.anglein for geom in self.xrfgeometries]

    @xrf_anglesin.setter
    def xrf_anglesin(self, value):
        for geom, pos in zip(self.xrfgeometries, value):
            geom.anglein = pos

    @property
    def xrf_anglesout(self):
        return [geom.angleout for geom in self.xrfgeometries]

    @xrf_anglesout.setter
    def xrf_anglesout(self, value):
        for geom, pos in zip(self.xrfgeometries, value):
            geom.angleout = pos

    @property
    def reference(self):
        return self._reference

    @reference.setter
    def reference(self, value):
        value = units.Quantity(value)
        try:
            value = value.to("Hz")
        except pinterrors.DimensionalityError:
            value = value.to("dimensionless")
        except pinterrors.DimensionalityError:
            raise ValueError(
                "Reference should be in Hz (flux) or dimensionless (counts)"
            )
        self._reference = value

    @property
    def defaultexpotime(self):
        return self._defaultexpotime

    @defaultexpotime.setter
    def defaultexpotime(self, value):
        if value is None:
            value = np.nan
        self._defaultexpotime = units.Quantity(value, "s")

    @property
    def defaultreferencetime(self):
        return self._defaultreferencetime

    @defaultreferencetime.setter
    def defaultreferencetime(self, value):
        if value is None:
            self._defaultreferencetime = None  # take the time from the scan
        else:
            self._defaultreferencetime = units.Quantity(value, "s")

    def fluxtocps(self, energy, flux, weights=None):
        return self.diodeI0.fluxtocps(energy, flux, weights=weights).magnitude

    def responsetoflux(self, energy, response, weights=None, time=None):
        return self.diodeI0.responsetoflux(
            energy, response, weights=weights, time=time
        ).magnitude

    def xrfnormop(
        self, energy, expotime=None, reference=None, referencetime=None, weights=None
    ):
        """
        Returns a function to be applied to I0 before normalizing XRF signal.

        Args:
            energy(num|array): source lines (keV)
            expotime(Optional(num)): original exposure time (sec)
            reference(Optional(num)): iodet (counts) or flux (photons/sec)
                                      to which the data should be normalized
            referencetime(Optional(num)): time to which the data should be normalized
            weights(Optional(num|array)): source line weights

        Returns:
            op(linop): raw diode conversion operator
            Fref(num): flux in photons/s to which the data
                       is normalized after data/op(diode)
            tref(num): time in s to which the data is normalized
                       after data/op(diode)
            t(num):    time in s of raw data
        """
        if expotime is None:
            expotime = self.defaultexpotime
        else:
            expotime = units.Quantity(expotime, "s")
        if reference is None:
            reference = self.reference
        else:
            reference = units.Quantity(reference, "dimensionless")
        if referencetime is None:
            referencetime = self.defaultreferencetime
        else:
            referencetime = units.Quantity(referencetime, "s")
        ret = self.diodeI0.xrfnormop(
            energy, expotime, reference, referencetime=referencetime, weights=weights
        )
        return ret + (units.umagnitude(expotime, "s"),)

    def quantinfo(self, *args, **kwargs):
        """
        Args:
            see xrfnormop

        Returns:
            op(linop): raw diode conversion operator
            info(list(dict)): pymca quantification info (flux and geometry)
        """
        xrfnormop, flux, time, expotime = self.xrfnormop(*args, **kwargs)

        quants = []
        for geom in self.xrfgeometries:
            quant = {"flux": flux, "time": time}
            quant["activearea"] = geom.detector.activearea
            quant["anglein"] = geom.anglein
            quant["angleout"] = geom.angleout
            quant["distance"] = self.getxrfdistance()
            quants.append(quant)

        return xrfnormop, quants

    def _beamfilter_transmission(self, energy, weights=None):
        Tbeamfilters = [
            weightedsum(
                geom.detector.filter_transmission(energy, source=True), weights=weights
            )
            for geom in self.xrfgeometries
        ]
        if len(set(Tbeamfilters)) > 1:
            beamfilters = "\n".join(
                [
                    " Geometry {}: {}".format(i, geom.beamfilters())
                    for i, geom in enumerate(self.xrfgeometries, 1)
                ]
            )
            raise RuntimeError(
                "Geometries should have the same beam filters:\n{}".format(beamfilters)
            )
        else:
            Tbeamfilters = Tbeamfilters[0]
        return Tbeamfilters

    def addsamplecovers(self, attenuators):
        for geom in self.xrfgeometries:
            for att in attenuators:
                geom.detector.addsamplecover(*att)

    def adddetectorfilters(self, attenuators):
        for geom in self.xrfgeometries:
            for att in attenuators:
                geom.detector.adddetectorfilter(*att)

    def addbeamfilters(self, attenuators):
        for geom in self.xrfgeometries:
            for att in attenuators:
                geom.detector.addbeamfilter(*att)

    def addtransmissionfilters(self, attenuators):
        for att in attenuators:
            name, material, thickness = att
            # This should be a beam filter so ignore name
            self.diodeIt.addbeamfilter(None, material, thickness)

    def I0op(self, energy, expotime=None, weights=None, removebeamfilters=False):
        """
        Calculate the flux before the sample from the I0 diode response

        Args:
            energy(num|array): keV
            expotime(Optional(num)): sec
            weights(Optional(num|array)): source lines ratio's
            removebeamfilters(Optional(bool)): flux before beam filter
                                               instead of flux before sample
        """
        if expotime is None:
            expotime = self.defaultexpotime
        op = self.diodeI0.fluxop(energy, 0, time=expotime, weights=weights)

        if removebeamfilters:
            Tbeamfilters = self._beamfilter_transmission(energy, weights=weights)
            op = LinearOperator(1.0 / Tbeamfilters, 0) * op

        op.m = op.m.to("hertz")
        op.b = op.b.to("hertz")
        return op, expotime

    def Itop(self, energy, expotime=None, weights=None, removebeamfilters=False):
        """
        Calculate the flux after the sample from the It diode response

        Args:
            energy(num|array): keV
            expotime(Optional(num)): sec
            weights(Optional(num|array)): source lines ratio's
            removebeamfilters(Optional(bool)): flux before beam filter
                                               instead of flux before sample
        """
        if expotime is None:
            expotime = self.defaultexpotime
        op = self.diodeIt.fluxop(energy, 0, time=expotime, weights=weights)

        if removebeamfilters:
            # Filters after the sample (like atmosphere) are already taken into account
            Tbeamfilters = self._beamfilter_transmission(energy, weights=weights)
            op = LinearOperator(1.0 / Tbeamfilters, 0) * op

        op.m = op.m.to("hertz")
        op.b = op.b.to("hertz")
        return op, expotime

    def batchcalibrate_diodes(self, params_fixed, params_var):
        """
        Execute `calibrate_diodes` from multiple parameter sets

        Args:
            params_fixed(dict): common for each set
            params_var(list(dict)): overrides common for each set
        """
        for params_var_ in params_var:
            params = params_fixed.copy()
            params.update(params_var_)
            self.calibrate_diodes(**params)

    def _calibrate_prepare(self, params):
        # Get data
        base = params.pop("base", None)
        sample = params.pop("sample", None)
        dataset = params.pop("dataset", None)
        specnr = params.pop("specnr", None)
        if specnr is not None:
            if sample:
                sampledataset = "{}_{}".format(sample, dataset)
                specfile = os.path.join(
                    base, sample, sampledataset, "{}.dat".format(sampledataset)
                )
            else:
                specfile = os.path.join(base, dataset)

            # Data in spec file
            data, staticdata = self._parse_spec(specfile, specnr)
        else:
            # Data in params
            data, staticdata = self._parse_default(params)

        # Prepare devices
        resetdevices = params.pop("resetdevices", False)
        params2 = self._validate_data(
            data, staticdata, resetdevices=resetdevices, dark=params["dark"]
        )
        gaindiodeI0 = params.pop("gaindiodeI0", None)
        gaindiodeIt = params.pop("gaindiodeIt", None)
        self._set_diode_gain(gaindiodeI0, params2.get("gaindiodeI0", None), "diodeI0")
        self._set_diode_gain(gaindiodeIt, params2.get("gaindiodeIt", None), "diodeIt")

        return data

    def _calibrate_dark(self, params, data):
        docalibrate = params.pop("calibrate", True)
        fixdark = params.pop("fixdark", False)
        params.pop("fluxmin", None)
        params.pop("fluxmax", None)
        if docalibrate and not fixdark:
            for k, attr in zip(self._calibrate_fields_id(), ["diodeI0", "diodeIt"]):
                deviceinstance = getattr(self, attr)
                kcurrent = "{}_current".format(k)
                kcps = "{}_cps".format(k)
                if kcurrent in data and kcps in data and False:
                    deviceinstance.calibratedark(data[kcurrent])
                    deviceinstance.calibrateF0(data[kcps])
                elif kcps in data:
                    deviceinstance.calibratedark(data[kcps])

    def _calibrate_diode_gain(self, params, data):
        fitinfo = {}
        docalibrate = params.pop("calibrate", True)
        if docalibrate:
            energy = data["energy"].to("keV").magnitude
            for name in ["current", "cps"]:
                I0, It = self._calibrate_fields_name(name)
                if It in data and I0 in data:
                    flux = self.diodeIt.responsetoflux(energy, data[It])
                    fitinfo = self.diodeI0.calibrate(data[I0], flux, energy, **params)
                    break
        return fitinfo

    def calibrate_diodes(self, **params):
        """
        Calibrate diodes at a particular energy
        """
        # Get data and prepare devices
        data = self._calibrate_prepare(params)

        # Calibrate
        dark = params.pop("dark", False)
        plot = params.pop("plot", False)
        oscillator = params.pop("oscillator", False)
        if oscillator:
            if plot:
                self._show_oscillator_calib(data)
        else:
            if dark:
                self._calibrate_dark(params, data)
                if plot:
                    self._show_dark_calib(data)
            else:
                fitinfo = self._calibrate_diode_gain(params, data)
                if plot:
                    fluxmin = params.get("fluxmin", 0)
                    fluxmax = params.get("fluxmax", np.inf)
                    self._show_flux_calib(
                        data, fluxmin=fluxmin, fluxmax=fluxmax, fitinfo=fitinfo
                    )
        if plot:
            plt.show()

    def _set_diode_gain(self, gainasked, gaindata, attr):
        deviceinstance = getattr(self, attr)
        if gainasked is None:
            if gaindata is None:
                raise RuntimeError(
                    "Gain of diode {} must be specified".format(
                        deviceinstance.__class__.__name__
                    )
                )
            deviceinstance.gain = gaindata
        else:
            if gaindata is not None:
                gainasked = units.quantity_like(gainasked, gaindata)
                if not np.isclose(gainasked, gaindata):
                    raise RuntimeError(
                        "Data was collected with diode {} gain {}, but {} was expected".format(
                            gaindata, deviceinstance.__class__.__name__, gainasked
                        )
                    )
            deviceinstance.gain = gainasked

    def _check_device(self, attr, clsfactory, clsnameused, resetdevice=False):
        if clsnameused is None:
            return  # usage is not specified

        deviceinstance = getattr(self, attr)
        if clsnameused:
            # Device was used
            if deviceinstance is None:
                # Device was not expected to be used
                if resetdevice:
                    setattr(self, attr, clsnameused)
                else:
                    raise RuntimeError(
                        "Device {} was not expected to be used.".format(clsnameused)
                    )
            else:
                # Device was expected to be used
                deviceclass = clsfactory(clsnameused)
                if not isinstance(deviceinstance, deviceclass):
                    raise RuntimeError(
                        "Data was collected with diode {}, while diode {} was expected.".format(
                            deviceclass.__name__, deviceinstance.__class__.__name__
                        )
                    )
        elif clsnameused is None:
            # No idea whether device was used or not
            pass
        else:
            # Device was not used
            if deviceinstance is not None:
                # Device was expected to be used
                if resetdevice:
                    setattr(self, attr, clsnameused)
                else:
                    raise RuntimeError(
                        "Data was not collected with device {}.".format(
                            deviceinstance.__class__.__name__
                        )
                    )

    def _calibrate_fields_id(self):
        return [self.fluxid, self.transmissionid]

    def _calibrate_fields_name(self, name):
        return ["{}_{}".format(k, name) for k in self._calibrate_fields_id()]

    def _calibrate_fields(self):
        return (
            self._calibrate_fields_name("counts")
            + self._calibrate_fields_name("cps")
            + self._calibrate_fields_name("current")
            + self._calibrate_fields_name("photons")
            + self._calibrate_fields_name("flux")
            + ["time", "energy"]
        )

    def _parse_spec(self, specfile, specnr):
        fspec = spec.spec(specfile)

        # Motors
        motornames = self.instrument.motornames()
        specmotornames = self.instrument.specmotornames(motornames)
        motorvalues = fspec.getmotorvalues(specnr, specmotornames)
        staticdata = dict(zip(motornames, motorvalues))

        # Counters
        fields = self._calibrate_fields()
        labels = [self.instrument.speccounternames.get(k, "") for k in fields]
        ascounter = fspec.haslabels(specnr, labels)
        asmotor = fspec.hasmotors(specnr, labels)

        labels1 = [label for label, v in zip(labels, ascounter) if v]
        ufields1 = [name for name, v in zip(fields, ascounter) if v]
        data = fspec.getdata2(specnr, labels1).T
        data = {
            k: units.Quantity(v, self.instrument.units[k])
            for k, v in zip(ufields1, data)
        }

        labels2 = [label for label, v in zip(labels, asmotor) if v]
        ufields2 = [name for name, v in zip(fields, asmotor) if v]
        motorvalues = fspec.getmotorvalues(specnr, labels2)

        data.update(
            {
                k: units.Quantity(v, self.instrument.units[k])
                for k, v in zip(ufields2, motorvalues)
            }
        )
        return data, staticdata

    def _parse_default(self, params):
        fields = self._calibrate_fields()
        data = {}
        for k in fields:
            if k in params:
                data[k] = units.Quantity(params.pop(k), self.instrument.units[k])
        return data, params.pop("motors", {})

    def _validate_data(self, data, staticdata, resetdevices=False, dark=False):
        fields = self._calibrate_fields()

        # Add units where missing and fill with static data if needed
        for k in fields:
            if k in data:
                data[k] = units.Quantity(data[k], self.instrument.units[k])
            elif k in staticdata:
                data[k] = units.Quantity(staticdata[k], self.instrument.units[k])

        # Fill in missing exposure time
        if "time" not in data:
            data["time"] = self.defaultexpotime
        data["time"] = data["time"].to("s")

        # Check units of energy
        if "energy" in data:
            data["energy"] = data["energy"].to("keV")

        # Diode response: ensure units
        for k in self._calibrate_fields_id():
            # counts are not used anywhere, just needed to calculate cps
            kcounts = "{}_counts".format(k)
            kphotons = "{}_photons".format(k)
            kcurrent = "{}_current".format(k)
            kcps = "{}_cps".format(k)
            kflux = "{}_flux".format(k)

            if kcounts in data:
                data[kcounts] = data[kcounts].to("dimensionless")
                if kcps not in data:
                    data[kcps] = data[kcounts].to("dimensionless") / data["time"]

            if kcps in data:
                data[kcps] = data[kcps].to("Hz")

            if kphotons in data:
                data[kphotons] = data[kphotons].to("dimensionless")
                if kflux not in data:
                    data[kflux] = data[kphotons].to("dimensionless") / data["time"]

            if kflux in data:
                data[kflux] = data[kflux].to("Hz")

            if kcurrent in data:
                data[kcurrent] = data[kcurrent].to("A")

        # Verify, create or replace devices
        I0, It = self._calibrate_fields_name("cps")
        diodeI0name = self.instrument.diodeI0(staticdata)
        diodeItname = self.instrument.diodeIt(staticdata)
        opticsname = self.instrument.optics(staticdata)
        self._check_device(
            "optics",
            xrayoptics.XrayOptics.clsfactory,
            opticsname,
            resetdevice=resetdevices,
        )
        self._check_device(
            "diodeI0", diodes.PNdiode.clsfactory, diodeI0name, resetdevice=resetdevices
        )
        self._check_device(
            "diodeIt", diodes.PNdiode.clsfactory, diodeItname, resetdevice=resetdevices
        )

        # Diode gains from response
        ret = {}
        if "energy" in data and not dark:
            energy = data["energy"].to("keV").magnitude
            for k, attr in zip(self._calibrate_fields_id(), ["diodeI0", "diodeIt"]):
                kcps = "{}_cps".format(k)
                if kcps in data:
                    for name in ["current", "flux"]:
                        j = "{}_{}".format(k, name)
                        if j in data:
                            ret["gain" + attr] = getattr(self, attr).gainfromresponse(
                                data[kcps], data[j], energy=energy
                            )
                            break
        return ret

    def _show_dark_calib(self, data):
        f, axs = plt.subplots(1, 2)
        color = None
        for k, attr, ax in zip(
            self._calibrate_fields_name("cps"), ["diodeI0", "diodeIt"], axs
        ):
            self._plot_dark_calib(attr, data[k], attr, ax, color=color)
            color = next(ax._get_lines.prop_cycler)["color"]
        plt.tight_layout()

    def _plot_dark_calib(self, deviceattr, response, name, ax, color=None):
        deviceinstance = getattr(self, deviceattr)
        ax.plot(response, "o", label="spec", color=color)
        ax.axhline(
            y=deviceinstance.fluxtocps(5, 0).magnitude, label="calc", color=color
        )
        ax.set_xlabel("points")
        ax.set_ylabel("{} ({:~})".format(name, response.units))
        ax.set_title("{:~.1e}".format(deviceinstance.gain))
        ax.legend(loc="best")

    def _plot_diode_flux(self, attr, energy, response, ax, color=None):
        deviceinstance = getattr(self, attr)

        if response.size == 1:
            response = units.asqarray([response[0] * 0.9, response[0] * 1.1])

        flux = deviceinstance.responsetoflux(energy, response)

        lines = ax.plot(flux, response, label=attr, color=color)
        ax.set_xlabel("flux@sample (ph/s)")
        ax.set_ylabel("{} ({:~})".format(attr, response.units))
        ax.set_title("{:~.1e}".format(deviceinstance.gain))
        ax.legend(loc="best")

        return lines[0].get_color()

    def _show_flux_calib(self, data, fluxmin=0, fluxmax=np.inf, fitinfo={}):
        f, (ax1, ax2) = plt.subplots(1, 2)

        I0, It = self._calibrate_fields_name("cps")
        _, fluxt = self._calibrate_fields_name("flux")

        energy = data["energy"].to("keV").magnitude
        It = units.asqarray(data[It])
        if fluxt in data:
            flux = units.asqarray(data[fluxt])
            ax1.plot(flux, It, "o", label="spec")
        color1 = self._plot_diode_flux("diodeIt", energy, It, ax1)
        color2 = next(ax1._get_lines.prop_cycler)["color"]

        xoff = 0.3
        yoff = 0.01
        info = self.diodeIt.fluxcpsinfo(energy)
        text = "\n".join(["{} = {}".format(k, v) for k, v in info.items()])
        ax1.text(xoff, yoff, text, transform=ax1.transAxes, verticalalignment="bottom")

        flux = self.diodeIt.responsetoflux(energy, It)
        fluxmin = units.Quantity(fluxmin, "hertz").to("hertz")
        fluxmax = units.Quantity(fluxmax, "hertz").to("hertz")
        I0 = units.asqarray(data[I0])
        indfit = (flux < fluxmin) | (flux > fluxmax)
        if any(indfit):
            ax2.plot(flux[indfit], I0[indfit], "o", color="#cccccc")
        indfit = np.logical_not(indfit)
        if any(indfit):
            ax2.plot(flux[indfit], I0[indfit], "o", label="diodeIt", color=color1)

        self._plot_diode_flux("diodeI0", energy, I0, ax2, color=color2)

        if fitinfo:
            text = "\n".join(["{} = {}".format(k, v) for k, v in fitinfo.items()])
            ax2.text(
                xoff, yoff, text, transform=ax2.transAxes, verticalalignment="bottom"
            )

        plt.tight_layout()

    def _plot_oscillator(self, attr, current, ax, color=None):
        deviceinstance = getattr(self, attr)

        op = deviceinstance.op_currenttocps()
        cps = op(current).to("Hz")
        current = current.to("pA")

        lines = ax.plot(current, cps, label=attr, color=color)
        ax.set_xlabel("{} ({:~})".format(attr, current.units))
        ax.set_ylabel("{} ({:~})".format(attr, cps.units))
        ax.set_title("{:~.1e}, {:~.1e}".format(op.m.to("Hz/pA"), op.b.to("Hz")))
        ax.legend(loc="best")

        return lines[0].get_color()

    def _show_oscillator_calib(self, data):
        f, (ax1, ax2) = plt.subplots(1, 2)

        I0, It = self._calibrate_fields_name("cps")
        cur0, curt = self._calibrate_fields_name("current")
        I0, It = data[I0], data[It]
        cur0, curt = data[cur0], data[curt]

        lines = ax1.plot(curt.to("pA"), It.to("Hz"), "o", label="spec")
        color1 = lines[0].get_color()
        color2 = next(ax1._get_lines.prop_cycler)["color"]
        self._plot_oscillator("diodeIt", curt, ax1, color=color1)

        ax2.plot(cur0.to("pA"), I0.to("Hz"), "o", label="spec", color=color2)
        self._plot_oscillator("diodeI0", cur0, ax2, color=color2)


def list_detectors():
    return list(xrfdetectors.registry.keys())


def list_sources():
    return list(xraysources.registry.keys())


def list_optics():
    return list(xrayoptics.registry.keys())


def list_geometries():
    return list(xrfgeometries.registry.keys())


def list_diodes():
    return list(diodes.registry.keys())


def list_instruments():
    return list(configuration.registry.keys())


def print_devices():
    print("instrument = {}".format(list_instruments()))
    print("xrfgeometry = {}".format(list_geometries()))
    print("xrfdetector = {}".format(list_detectors()))
    print("diodeI0,diodeIt = {}".format(list_diodes()))
    print("source = {}".format(list_sources()))
    print("optics = {}".format(list_optics()))


class SXM(QXRFGeometry):
    def __init__(self, **kwargs):
        kwargs["diodeIt"] = kwargs.get("diodeIt", "idet")
        kwargs["xrfdetector"] = kwargs.get("xrfdetector", "leia")
        kwargs["xrfgeometry"] = kwargs.get("xrfgeometry", "sxm120")
        kwargs["instrument"] = kwargs.get("instrument", "sxm")
        kwargs["simplecalibration"] = kwargs.get("simplecalibration", True)
        super(SXM, self).__init__(**kwargs)


class SXM1(QXRFGeometry):
    aliases = ["sxm1_kb", "sxm1_kb"]

    def __init__(self, **kwargs):
        kwargs["diodeI0"] = kwargs.get("diodeI0", "iodet1")
        kwargs["diodeIt"] = kwargs.get("diodeIt", "idet")
        kwargs["optics"] = kwargs.get("optics", "kb")
        kwargs["xrfdetector"] = kwargs.get("xrfdetector", "leia")
        kwargs["xrfgeometry"] = kwargs.get("xrfgeometry", "sxm120")
        kwargs["instrument"] = kwargs.get("instrument", "sxm")
        kwargs["simplecalibration"] = kwargs.get("simplecalibration", True)
        super(SXM1, self).__init__(**kwargs)


class SXM2(QXRFGeometry):
    aliases = ["sxm2_kb", "sxm2_kb"]

    def __init__(self, **kwargs):
        kwargs["diodeI0"] = kwargs.get("diodeI0", "iodet2")
        kwargs["diodeIt"] = kwargs.get("diodeIt", "idet")
        kwargs["optics"] = kwargs.get("optics", "kb")
        kwargs["xrfdetector"] = kwargs.get("xrfdetector", "leia")
        kwargs["xrfgeometry"] = kwargs.get("xrfgeometry", "sxm120")
        kwargs["instrument"] = kwargs.get("instrument", "sxm")
        kwargs["simplecalibration"] = kwargs.get("simplecalibration", True)
        super(SXM2, self).__init__(**kwargs)


class SXM1_Unfocused(QXRFGeometry):
    def __init__(self, **kwargs):
        kwargs["diodeI0"] = kwargs.get("diodeI0", "iodet1")
        kwargs["diodeIt"] = kwargs.get("diodeIt", "idet")
        kwargs["xrfdetector"] = kwargs.get("xrfdetector", "leia")
        kwargs["xrfgeometry"] = kwargs.get("xrfgeometry", "sxm120")
        kwargs["instrument"] = kwargs.get("instrument", "sxm")
        kwargs["simplecalibration"] = kwargs.get("simplecalibration", True)
        super(SXM1_Unfocused, self).__init__(**kwargs)


class SXM2_Unfocused(QXRFGeometry):
    def __init__(self, **kwargs):
        kwargs["diodeI0"] = kwargs.get("diodeI0", "iodet2")
        kwargs["diodeIt"] = kwargs.get("diodeIt", "idet")
        kwargs["xrfdetector"] = kwargs.get("xrfdetector", "leia")
        kwargs["xrfgeometry"] = kwargs.get("xrfgeometry", "sxm120")
        kwargs["instrument"] = kwargs.get("instrument", "sxm")
        kwargs["simplecalibration"] = kwargs.get("simplecalibration", True)
        super(SXM2_Unfocused, self).__init__(**kwargs)


factory = QXRFGeometry.factory
registry = QXRFGeometry.clsregistry

# -*- coding: utf-8 -*-

import math

# import warnings
import logging
import collections

from ..patch.pint import ureg
from ..math import linop
from ..math import fit1d
from ..utils import units
from ..materials import compoundfromformula
from ..materials import compoundfromname
from ..materials import multilayer
from ..materials import element
from ..resources import resource_filename
from ..utils import constants
from ..math.utils import round_sig
from ..geometries import diode as diodegeometries
from ..optics import xray as xrayoptics
from ..sources import xray as xraysources
from ..simulation.classfactory import with_metaclass
from ..math import noisepropagation
from ..utils import instance
from . import base
from ..utils import lut
from ..materials.utils import reshape_spectrum_lines

import numpy as np
import silx.math.fit as fit
import scipy.interpolate
import matplotlib.pyplot as plt
import scipy.optimize
from pint import errors as pinterrors
from collections import Counter

logger = logging.getLogger(__name__)

#
#
# Unit of elementary charge: 1 e = 1.6e-19 C
# Electronvolt: 1 eV = 1.6e-19 J
# 1 A/W = 1 C/J = 1.6e-19 C/eV = 1 e/eV
#
#
# Oscillator: current (A) to counts/sec (Hz)
#  I(Hz) = Fmax(Hz)/Vmax(V) . Vmax(V)/Imax(A) . I(A) + F0(Hz)
#    SXM: Vmax/Imax set
#    MICRODIFF: Imax set
#
# PNdiode: flux (ph/s) to current(A)
#  I(A) = I(ph/s).E(eV/ph).1(e).(1-T)/Ehole(eV) + D(A)
#       = P(W).1(e).(1-T)/Ehole(eV) + D(A)
#       = P(W).R(e/eV) + D(A)
#    T: transmission of the diode
#    P: radiative power
#    e: elementary charge
#    R: spectral responsivity
#
# Sample flux:
#  Is(ph/s) = I0(ph/s).Ts
#  I(ph/s) = I0(ph/s).Y = Is(ph/s).Y/Ts
#    I: flux seen by the diode
#    I0: the flux of the source
#    Is: flux seen by the sample
#    Y: yield of secondary target (includes solid angle of detector but not attenuation)
#    Ts: transmission of the source lines (optics,diode/target,filters)
#
# For multiple source lines (i) and multiple secondary lines (j):
#   Iis(ph/s) = I0(ph/s).wi.Tis
#   Is(ph/s) = sum_i[Iis] = I0(ph/s) . sum_i [wi.Tis]
#   Iij(ph/s) = I0(ph/s).wi.Yij = Is(ph/s).wi.Yij/sum_i[wi.Tis]
#
#   I(A) = SUM_ij [Iij(ph/s).Ej(eV/ph).1(e).(1-Tj)/Ehole(eV)] + D(A)
#        = SUM_ij [Iij(ph/s).Cj] + D(A)
#        = Is(ph/s) . SUM_i [SUM_j[wi.Yij.Cj]] / SUM_k [wk.Tks] + D(A)
#        = Is(ph/s) . Cs + D(A)
#
#   wi: source line fraction
#   Cj (C/ph): charge per photon hitting the diode
#   Cs (C/ph): charge per photon hitting the sample
#
# The source fraction at the sample position:
#   wis = Iis/sum_i[Iis] = wi.Tis / sum_i[wi.Tis]
#
# An absolute diode measures both I(A) and P(W) so that
# another diode which only measures I(A) can be calibrated:
#  I_abs(A) = P(W).R_abs(e/eV) + D_abs(A)
#  I(A) = P(W).FACT.R_abs(e/eV) + D(A)
#
#   FACT = R(e/eV)/R_abs(e/eV)
#        = [I(A)-D(A)]/[I_abs(A)-D_abs(A)]
#        = (1-T)/Ehole(eV) . Ehole_abs(eV)/(1-T_abs)
#


class GainRounder(object):
    def __init__(self, m=None, base=None):
        self.m = m
        self.base = base

    def __getstate__(self):
        return {"m": self.m, "base": self.base}

    def __setstate__(self, state):
        self.m = state["m"]
        self.base = state["base"]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.m == other.m and self.base == other.base
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __call__(self, x):
        if self.m is None and self.base is None:
            return x
        gain = x.magnitude
        u = x.units
        if self.m is not None:
            gain = gain / self.m
        if self.base is not None:
            gain = math.log(gain, self.base)
        gain = int(round(gain))
        if self.base is not None:
            gain = self.base**gain
        if self.m is not None:
            gain = self.m * gain
        return units.Quantity(gain, u)


class Oscillator(object):
    """This class describes a voltage-to-frequency convertor"""

    def __init__(self, Fmax=None, F0=None, Vmax=None):
        """
        Args:
            Fmax(num): maximal frequency (Hz)
            F0(num): frequency offset (Hz)
            Vmax(num): voltage corresponding to Fmax (V)
        """
        self.Fmax = Fmax
        self.F0 = F0
        self.Vmax = Vmax

    def __getstate__(self):
        return {"Fmax": self.Fmax, "F0": self.F0, "Vmax": self.Vmax}

    def __setstate__(self, state):
        self.Fmax = state["Fmax"]
        self.F0 = state["F0"]
        self.Vmax = state["Vmax"]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (
                self.Fmax == other.Fmax
                and self.F0 == other.F0
                and self.Vmax == other.Vmax
            )
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def Fmax(self):
        return self._Fmax

    @Fmax.setter
    def Fmax(self, value):
        value = units.Quantity(value, "hertz")
        self._Fmax = value.to("hertz")

    @property
    def F0(self):
        return self._F0

    @F0.setter
    def F0(self, value):
        value = units.Quantity(value, "hertz")
        self._F0 = value.to("hertz")

    @property
    def Vmax(self):
        return self._Vmax

    @Vmax.setter
    def Vmax(self, value):
        value = units.Quantity(value, "volt")
        self._Vmax = value.to("volt")

    def tojson(self):
        kwargs = {"Fmax": self.Fmax, "F0": self.F0, "Vmax": self.Vmax}
        return {"classname": self.__class__.__name__, "kwargs": kwargs}

    def __str__(self):
        return "y hertz = {:~e}/{:~} * x V + {:~}".format(self.Fmax, self.Vmax, self.F0)

    def op_voltagetocps(self):
        """Operator to convert voltage to counts-per-second

        Args:
            None

        Returns:
            callable: slope (cps/V), intercept(cps)
        """
        return linop.LinearOperator(self.Fmax / self.Vmax, self.F0)

    def op_cpstovoltage(self):
        """Operator to convert counts-per-second to voltage

        Args:
            None

        Returns:
            callable: slope (V/cps), intercept(V)
        """
        return self.op_voltagetocps().inverse


class PNdiode(with_metaclass(base.SolidState)):

    ELCHARGE = ureg.Quantity(1, ureg.elementary_charge)
    # ureg.Quantity(1.6e-19,ureg.coulomb) # approximation used in spec

    def __init__(
        self,
        gain=None,
        gainrounder=None,
        darkcurrent=None,
        oscillator=None,
        secondarytarget=None,
        simplecalibration=False,
        optics=None,
        Vmax=None,
        beforesample=None,
        **kwargs
    ):
        """
        Args:
            gain(num): V/A (default) or A
            gainrounder(callable): discrete gain settings
            darkcurrent(num): C/s
            oscillator(Oscillator):
            secondarytarget(multilayer): optional secondary target
            simplecalibration(bool): simplify
            optics:
            Vmax: when oscillator is missing
        """
        if Vmax is None and oscillator is not None:
            Vmax = oscillator.Vmax
        self.Vmax = Vmax
        self.oscillator = oscillator
        self.gainrounder = gainrounder
        self.gain = gain
        self.darkcurrent = darkcurrent
        self.beforesample = beforesample
        self.optics = optics
        self.secondarytarget = secondarytarget
        self.simplecalibration = simplecalibration
        self._lut_chargepersamplephoton = lut.LUT()
        super(PNdiode, self).__init__(**kwargs)

    def __getstate__(self):
        state = super(PNdiode, self).__getstate__()
        state["oscillator"] = self.oscillator
        state["Vmax"] = self.Vmax
        state["gainrounder"] = self.gainrounder
        state["gain"] = self.gain
        state["darkcurrent"] = self.darkcurrent
        state["beforesample"] = self.beforesample
        state["optics"] = self.optics
        if self.secondarytarget is not None:
            state["secondarytarget"] = self.secondarytarget
        state["simplecalibration"] = self.simplecalibration
        state["_lut_chargepersamplephoton"] = self._lut_chargepersamplephoton
        return state

    def __setstate__(self, state):
        super(PNdiode, self).__setstate__(state)
        self.oscillator = state["oscillator"]
        self.Vmax = state["Vmax"]
        self.gainrounder = state["gainrounder"]
        self.gain = state["gain"]
        self.darkcurrent = state["darkcurrent"]
        self.beforesample = state["beforesample"]
        self.optics = state["optics"]
        self.secondarytarget = state.get("secondarytarget", None)
        self.simplecalibration = state["simplecalibration"]
        self._lut_chargepersamplephoton = state["_lut_chargepersamplephoton"]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if not (super(PNdiode, self).__eq__(other)):
                return False
            if not self._eq_secondarytarget(other):
                return False
            return (
                self.oscillator == other.oscillator
                and self.Vmax == other.Vmax
                and self.gainrounder == other.gainrounder
                and self.gain == other.gain
                and self.darkcurrent == other.darkcurrent
                and self.beforesample == other.beforesample
                and self.optics == other.optics
                and self.simplecalibration == other.simplecalibration
                and self._lut_chargepersamplephoton == other._lut_chargepersamplephoton
            )
        else:
            return False

    def _eq_secondarytarget(self, other):
        if (self.secondarytarget is None) ^ (other.secondarytarget is None):
            return False
        if self.secondarytarget is not None:
            # Circular reference:
            gself = self.secondarytarget.geometry.detector
            gother = other.secondarytarget.geometry.detector
            self.secondarytarget.geometry.detector = None
            other.secondarytarget.geometry.detector = None
            eqtarget = self.secondarytarget == other.secondarytarget
            self.secondarytarget.geometry.detector = gself
            other.secondarytarget.geometry.detector = gother
            return eqtarget
        return True

    @property
    def optics(self):
        return self._optics

    @optics.setter
    def optics(self, optics):
        if optics:
            if not instance.isarray(optics):
                optics = [optics]

            def proc(dev):
                if instance.isstring(dev):
                    return xrayoptics.factory(dev)
                else:
                    return dev

            optics = [proc(dev) for dev in optics]
        else:
            optics = []
        self._optics = optics

    @property
    def caliboptic(self):
        for dev in self.optics:
            if dev.uselut:
                return dev
        raise RuntimeError("No optics with calibratable transmission (with LUT)")

    @property
    def simplecalibration(self):
        return self._simplecalibration and not self._lut_chargepersamplephoton.isempty()

    @simplecalibration.setter
    def simplecalibration(self, value):
        self._simplecalibration = value

    @property
    def gainrounder(self):
        return self._gainrounder

    @gainrounder.setter
    def gainrounder(self, value):
        if value is None:
            value = GainRounder()
        self._gainrounder = value

    @property
    def Rout(self):
        """Output resistance of the ammeter (V/A)"""
        try:
            return self.gain.to("V/A")
        except pinterrors.DimensionalityError:
            return self.Vmax / self.Imax

    @Rout.setter
    def Rout(self, value):
        value = units.Quantity(value, "V/A")
        try:
            self.gain = value.to(self.gainunits)
        except pinterrors.DimensionalityError:
            self.gain = (self.Vmax / value).to(self.gainunits)

    @property
    def Imax(self):
        try:
            return self.gain.to("A")
        except pinterrors.DimensionalityError:
            return self.Vmax / self.gain.to("V/A")

    @Imax.setter
    def Imax(self, value):
        value = units.Quantity(value, "A")
        try:
            self.gain = value.to(self.gainunits)
        except pinterrors.DimensionalityError:
            self.gain = (self.Vmax / value).to(self.gainunits)

    @property
    def gain(self):
        """Vmax/Imax or Imax"""
        return self._gain

    @gain.setter
    def gain(self, value):
        try:
            u = self.gainunits
        except AttributeError:
            try:
                gain = value.to("A")
            except pinterrors.DimensionalityError:
                gain = value.to("V/A")
        else:
            gain = units.Quantity(value, u).to(u)
        self._gain = self.gainrounder(gain)

    @property
    def gainunits(self):
        return self._gain.units

    @property
    def darkcurrent(self):
        return self._darkcurrent

    @darkcurrent.setter
    def darkcurrent(self, value):
        self._darkcurrent = units.Quantity(value, "A")

    def __str__(self):
        fmt = "PN-diode:\n{}\n"
        args = [super(PNdiode, self).__str__()]
        if self.simplecalibration:
            fmt += "LUT (e-/sample photon):\n {}\n"
            lut = self._lut_chargepersamplephoton.zip("keV", "e")
            lut = "\n ".join("{:~}: {:~}".format(*xy) for xy in lut)
            if not lut:
                lut = None
            args.append(lut)
        else:
            fmt += "Secondary target:\n {}\n"
            target = str(self.secondarytarget)
            target = "\n ".join(target.split("\n"))
            args.append(target)
            fmt += "Optics:\n {}\n"
            optics = "\n".join(str(dev) for dev in self.optics)
            optics = "\n ".join(optics.split("\n"))
            if not optics:
                optics = None
            args.append(optics)
            fmt += "Before sample: {}\n"
            args.append(self.beforesample)
        fmt += "Ammeter:\n Gain = {:~e}\n "
        args.append(self.gain)
        fmt += "Dark current = {:~}\n"
        args.append(self.darkcurrent.to("e/s"))
        if self.Vmax is not None:
            fmt += " Output voltage = {:~}\n"
            args.append(self.Vmax)
        fmt += "Voltage-to-Frequency:\n {}"
        args.append(self.oscillator)
        return fmt.format(*args)

    def _diode_absorbance(self, energy, thickness):
        # list self.absorbance but with customn thickness
        return self.material.mass_att_coeff(energy) * self.material.density * thickness

    def _diode_transmission(self, energy, thickness):
        # list self.transmission but with customn thickness
        return np.exp(-self._diode_absorbance(energy, thickness))

    def _diode_attenuation(self, energy, thickness):
        # list self.absorbance but with customn thickness
        return 1 - self._diode_transmission(energy, thickness)

    def _spectral_responsivity_infthick(self):
        """Generated current per radiative power for an infinitely thick diode
        Returns:
            num: A/W or e/eV
        """
        return (self.ELCHARGE / self.ehole.to("eV")).to("A/W")

    def spectral_responsivity(self, energy):
        """Generated current per radiative power

        Args:
            energy(num or array): keV

        Returns:
            num or array: A/W or e/eV
        """
        return self.attenuation(energy) * self._spectral_responsivity_infthick()

    def model_spectral_responsivity(self, energy, thickness, ehole):
        """Function to fit measured spectral responsivity"""
        return self._diode_attenuation(energy, thickness) / ehole

    def fit_spectral_responsivity(self, energy, response):
        """Calculate d and Ehole by fitting:

            I(A) = P(W) . 1(e) . (1-exp(-mu.rho.d))/Ehole(eV) + D(A)
            response(A/W) = (I(A)-D(A))/P(W) = 1(e).(1-exp(-mu.rho.d))/Ehole(eV)

        Args:
            energy(array): keV
            response(array): A/W

        Returns:
            None
        """
        # with warnings.catch_warnings():
        #    warnings.simplefilter("ignore")
        #    ehole = 1/np.max(response)
        #    muL = self.material.mass_att_coeff(energy)*self.material.density
        #    thickness = -np.log(1-(response*ehole).magnitude)/muL
        #    thickness = np.median(thickness[np.isfinite(thickness)])
        #    ehole = ehole.magnitude
        ehole = self.ehole.to("eV").magnitude
        thickness = self.thickness
        try:
            p, cov_matrix = fit.leastsq(
                self.model_spectral_responsivity, energy, response, [thickness, ehole]
            )
        except:
            logger.debug(
                "Could not fit spectral response of {}".format(self.__class__.__name__)
            )
            return
        self.thickness = p[0]  # "cm"
        self.ehole = units.Quantity(p[1], "eV")  # eV

    def plot_spectral_responsivity(self, energy):
        plt.plot(energy, self.spectral_responsivity(energy).to("A/W"))
        ax = plt.gca()
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Spectral responsivity (A/W)")
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        plt.title(self.__class__.__name__)

    def plot_response(self, energy, flux, weights=None, current=False):
        """
        Args:
            energy(num or array): source energies
                                  shape = [nSource x] nSourceLines
            flux(num or array): source energies
                                shape = [nSource x] nFluxes
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines
        """
        energy, weights, singlesource, singleline = reshape_spectrum_lines(
            energy, weights=weights
        )
        if current:
            response = (
                self.fluxtocurrent(energy, flux, weights=weights).to("A").magnitude
            )
            unit = "A"
        else:
            response = self.fluxtocps(energy, flux, weights=weights).to("Hz").magnitude
            unit = "Hz"
        for i, (x, y) in enumerate(zip(energy, response)):
            plt.plot(x, y, label="Source{}".format(i))
        # plt.axhline(y=self.oscillator.Fmax.to("Hz").magnitude,label="max")
        ax = plt.gca()
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Response ({})".format(unit))
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        plt.title(
            "{} @ {:.02e} ph/s, {:~.0e}".format(
                self.__class__.__name__, flux, self.gain
            )
        )
        if not singlesource:
            plt.legend()

    def _chargeperdiodephoton(self, energy):
        """
        Charge generated per photon hitting the diode: spectral responsivity multiplied by the energy

        Args:
            energy(num or array): keV

        Returns:
            num or array: C/ph
        """
        return self.spectral_responsivity(energy) * units.Quantity(energy, "keV")

    def _transmission_optics(self, energy):
        """
        Args:
            energy(array): source energies in keV (shape: nSource)

        Returns:
            array:
        """
        T = np.ones_like(energy)
        for dev in self.optics:
            T = T * dev.transmission(energy)
        return T

    def _source_transmission(self, energy):
        """
        Transmission between sample and point-before-detection (i.e. before source filter)

        Args:
            energy(array): keV

        Returns:
            array:
        """
        # Before sample:
        #  Direct detection  : point-before-detection - filter(source) - filter(det) - diode - optics - sample
        #  Indirect detection: point-before-detection - filter(source) - target - optics - sample
        #
        # After sample:
        #  Direct detection  : sample - optics - point-before-detection - filter(source) - filter(det) - diode
        #  Indirect detection: sample - optics - point-before-detection - filter(source) - target
        T = self._transmission_optics(energy)
        if self.beforesample:
            T = T * self.filter_transmission(energy, source=True)
            if self.secondarytarget is None:
                # Direct detection
                T = (
                    T
                    * self.transmission(energy)
                    * self.filter_transmission(energy, source=False)
                )
            else:
                # Indirect detection
                T = T * self.secondarytarget.transmission(energy)
        return T

    def _detection_rates(self, energy, weights):
        """
        Rate of line j at point-before-diode (i.e. without diode attenuation)
        multiplied by (and due to) rate of line i at point-before-detection (i.e. before source filter)

        Args:
            energy(array): source energies in keV
                           shape: [nSource x] nSourceLines
            weights(array): source rates
                            shape: [nSource x] nSourceLines

        Returns:
            energy2(array): energies (keV) of lines detected by the diode
                            shape: [nSource] x nSourceLines2
            rates(array): rates of lines detected by the diode
                          shape: [nSource x] nSourceLines2
        """
        energy, weights, singlesource, singleline = reshape_spectrum_lines(
            energy, weights=weights
        )
        if self.secondarytarget is None:
            # Directly detect the source spectrum
            wY = (
                weights
                * self.filter_transmission(energy, source=True)
                * self.filter_transmission(energy, source=False)
            )
            wY = wY.sum(axis=-1, keepdims=True)
        else:
            # Target spectrum energies and rates (ph/phsource)
            secondarytarget_ratedict = None
            if secondarytarget_ratedict is None:
                # includes source and detector filters
                # but not detector attenuation (withdetectorresponse=False)
                spectra = self.secondarytarget.xrayspectrum(
                    energy, weights=weights, withdetectorresponse=False, method="mc"
                )
                spectra = [dict(spectrum.spectrum()) for spectrum in spectra]
            else:
                # Soure weights before the sample:
                weights = weights * self.filter_transmission(energy, source=True)
                spectra = []
                for en, we in zip(energy, weights):
                    spectrum = Counter()
                    for e, w in zip(en, we):
                        rdict = secondarytarget_ratedict[e]
                        senergies = np.array(list(rdict.keys()))
                        # Product of:
                        #   source weight before the sample
                        #   rate after sample interaction
                        #   attenuation before detector
                        srates = (
                            w
                            * np.array(list(rdict.values()))
                            * self.filter_transmission(senergies, source=False)
                        )
                        for k, v in zip(senergies, srates):
                            spectrum[k] += v
                    spectra.append(dict(spectrum))
            # Spectra as arrays:
            energies = set()
            for spectrum in spectra:
                energies |= set(spectrum.keys())
            energy = []
            wY = []
            for spectrum in spectra:
                missing = energies - set(spectrum.keys())
                if missing:
                    spectrum.update(zip(missing, np.zeros(len(missing))))
                energy.append(list(spectrum.keys()))
                wY.append(list(spectrum.values()))
            energy = np.asarray(energy)
            wY = np.asarray(wY)
        if singlesource:
            return energy[0], wY[0]
        else:
            return energy, wY

    def _chargepersamplephoton(self, energy, weights=None, keepdims=False):
        """
        Charge generated per photon reaching the sample

        Args:
            energy(num or array): source energies in keV (shape: [nSource x] nSourceLines)
            weights(num or array): source line weights (shape: [nSource x] nSourceLines)

        Returns:
            num: C/ph (num or array(nSource [x 1]))
        """
        # Parse input (nSource x nSourceLines)
        energy, weights, singlesource, singleline = reshape_spectrum_lines(
            energy, weights=weights
        )
        kwargs = {"axis": -1, "keepdims": keepdims}
        if self.simplecalibration:
            # I(A)  = Is(ph/s) . Cs(C) + D(A)
            # Cs(C) = SUM_i [wi.LUTj]
            LUTj = self._lut_chargepersamplephoton(units.Quantity(energy, "keV"))
            # Sum over source lines:
            Cs = np.sum(LUTj * weights, **kwargs)  # nSource
        else:
            # I(A)  = Is(ph/s) . Cs(C) + D(A)
            # Cs(C) = SUM_i [SUM_j[wi.Yij.Cj]] / SUM_k [wk.Tks]
            #
            # wi, wk: source line fraction
            # Cj (C/ph): charge per photon hitting the diode
            # Cs (C/ph): charge per photon hitting the sample
            # Yij: rate of line j due to line i (including solid angle but not diode attenuation)
            # Tks: transmission of line k (point-before-detection <-> sample)
            Ts = self._source_transmission(energy)  # nSource x nSourceLines
            energy2, wiYij = self._detection_rates(
                energy, weights
            )  # nSource x nSecondary
            Cj = self._chargeperdiodephoton(energy2)  # nSource x nSecondary
            # Sum over secondary source lines:
            Cs = np.sum(wiYij * Cj, **kwargs) / np.sum(
                weights * Ts, **kwargs
            )  # nSource
        if singlesource:
            return Cs[0]
        else:
            return Cs

    def _calibrate_chargepersamplephoton(
        self,
        energy,
        Cscalib,
        weights=None,
        caliboption="optics",
        bound=False,
        plot=False,
    ):
        """
        Args:
            energy(num or array): source energies in keV (shape: nSourceLines)
            Cscalib(num): new charge-per-sample-photon (C or e)
            weights(num or array): source line weights (shape: nSourceLines)
            caliboption(str): "optics", "solidangle" or "thickness"

        Returns:
            Cscalc(num)
        """
        # Parse input (1 x nSourceLines)
        energy, weights, singlesource, singleline = reshape_spectrum_lines(
            energy, weights=weights
        )
        energy = energy[0]
        weights = weights[0]
        # I(A) = Is(ph/s).Cs + D(A)
        if self._simplecalibration:
            # Cs = SUM_i [SUM_j[wi.Yij.Cj]] / SUM_k [wk.Tks]
            #    = SUM_i [wi . Csi]
            self._lut_chargepersamplephoton.replace(
                units.Quantity(energy, "keV"), Cscalib / weights
            )
            Cscalc = self._chargepersamplephoton(energy, weights=weights)
        else:
            # Propagate correction to one of the components of Cs:
            #  Cs = SUM_i [SUM_j[wi.Yij.Cj]] / SUM_k [wk.Tks]
            #  Cj = (1-Tj).Cj_inf
            #  Cj_inf = Ej(eV/ph).1(e)/Ehole(eV)
            #  Yij = solidangle . ...
            #  Tj = 1-exp(-rho.mu(Ej).thickness)
            #
            #  Yij: rate of line j due to line i (including solid angle but not diode attenuation)
            #  Cj (C/ph): charge per photon hitting the diode
            #  Cs (C/ph): charge per photon hitting the sample
            if caliboption not in ["solidangle", "thickness", "optics"]:
                raise RuntimeError(
                    'caliboption must be "optics", "solidangle" or "thickness"'
                )
            if caliboption == "solidangle":
                Cs = self._chargepersamplephoton(energy, weights=weights)
                # solidangle' = solidangle*Cs'/Cs
                sa = self.geometry.solidangle * units.magnitude(
                    Cscalib / Cs, "dimensionless"
                )
                if sa <= 0 or sa > (2 * np.pi):
                    logger.warning(
                        "Diode solid angle of 4*pi*{} srad is not valid but will be accepted anyway".format(
                            sa / (4 * np.pi)
                        )
                    )
                self.geometry.solidangle = sa
                Cscalc = self._chargepersamplephoton(energy, weights=weights)
            else:
                # No analytical solution (in general): diode thickness or optics transmission
                # Fluorescence does not change (expensive to calculate)
                energy2, wiYij = self._detection_rates(
                    energy, weights
                )  # nSource x nSecondary, nSource x nSecondary
                # Solve: func(x) = (Cscalib-Cs(x))^2 = 0
                Cs, xorg, x0 = self._calibrate_init(
                    energy, caliboption=caliboption, bound=bound
                )
                # Remark: calling function Cs changes self based on first argument
                Cscalibm = Cscalib.to("e").magnitude

                def func(x):
                    return (
                        Cscalibm
                        - Cs(x, energy, weights, energy2, wiYij).to("e").magnitude
                    ) ** 2

                if bound:
                    x, info = scipy.optimize.bisect(
                        func, x0[0], x0[-1], full_output=True, disp=False
                    )
                    if not info.converged:
                        x = np.nan
                    if plot:
                        print(info)
                else:
                    x, infodict, ier, emsg = scipy.optimize.fsolve(
                        func, x0=x0, full_output=True
                    )
                    if plot:
                        print(infodict)
                        print(emsg, ier)
                    if ier != 1:
                        x = np.nan

                if plot:
                    plt.axvline(x0, linestyle="dashed")
                    plt.axvline(x)
                    plt.axhline(0)
                    xx = np.linspace(0.1, 1.5, 50)
                    yy = [func(xi) for xi in xx]
                    plt.plot(xx, yy)

                x = self._calibrate_apply(energy, x, xorg, caliboption=caliboption)
                Cscalc = Cs(x, energy, weights, energy2, wiYij)

                if plot:
                    print(bound, x0, x)
                    print(Cscalib, Cscalc)
                    plt.show()
        return Cscalc

    def _calibrate_init(self, energy, caliboption="optics", bound=False):
        if caliboption == "thickness":
            # Unknown x = diode thickness
            xorg = self.thickness
            if bound:
                attenuationmax = 0.999
                # use energy instead of energy2 because we need the maximal thickness
                thicknessmax = -np.log(1 - attenuationmax) / (
                    self.material.mass_att_coeff(energy) * self.material.density
                )
                x0 = (0.0, np.max(thicknessmax))
            else:
                x0 = xorg
            Cs = self._calibrate_thickness
        elif caliboption == "optics":
            # Unknown x = transmission
            xorg = self.caliboptic.transmission(energy)
            if len(energy) > 1:
                bound = False
            if bound:
                x0 = 0.0, 1.0
            else:
                x0 = xorg
            Cs = self._calibrate_optics
        else:
            raise RuntimeError('caliboption must be "optics" or "thickness"')
        return Cs, xorg, x0

    def _calibrate_apply(self, energy, x, xorg, caliboption="optics"):
        if np.isnan(x):
            logger.error("Diode calibration not succesfull")
            x = xorg
        if caliboption == "thickness":
            x = instance.asscalar(x)
            if x <= 0:
                logger.warning(
                    "Diode thickness of {} um is not valid but will be accepted anyway".format(
                        x * 1e-4
                    )
                )
            self.thickness = x
        elif caliboption == "optics":
            if x < 0 or x > 1:
                logger.warning(
                    "Transmission of {} % is not valid but will be accepted anyway".format(
                        x * 100
                    )
                )
            self.caliboptic.set_transmission(energy, x)
        else:
            raise RuntimeError('caliboption must be "optics" or "thickness"')
        return x

    def _calibrate_thickness(self, x, energy, weights, energy2, wiYij):
        self.thickness = instance.asscalar(x)
        Ts = self._source_transmission(energy)  # nSource
        Cj = self._chargeperdiodephoton(energy2)  # nSecondary
        return np.sum(wiYij * Cj[np.newaxis, :]) / np.sum(weights * Ts)

    def _calibrate_optics(self, x, energy, weights, energy2, wiYij):
        self.caliboptic.set_transmission(energy, x)
        Ts = self._source_transmission(energy)  # nSource
        Cj = self._chargeperdiodephoton(energy2)  # nSecondary
        return np.sum(wiYij * Cj[np.newaxis, :]) / np.sum(weights * Ts)

    def samplelineweights(self, energy, weights=None):
        """Source weights after transmission

        Args:
            energy(num or array): source energies in keV (shape: nSource)
            weights(num or array): source line weights (shape: nSource)

        Returns:
            weights(num or array): source line weights at the sample position
        """
        # Parse input
        energy, weights, singlesource, singleline = reshape_spectrum_lines(
            energy, weights=weights
        )
        if self.simplecalibration or not self.beforesample:
            return weights
        else:
            # wis = Iis/sum_i[Iis] = wi.Tis / sum_i[wi.Tis]
            weights *= self._source_transmission(energy)
            weights /= weights.sum(axis=-1, keepdims=True)
            return weights

    def op_currenttovoltage(self):
        """Operator to convert current to voltage

        Args:
            None

        Returns:
            callable: slope (V/A), intercept = 0V
        """
        return linop.LinearOperator(self.Rout, units.Quantity(0, "volt"))

    def op_voltagetocurrent(self):
        """Operator to convert voltage to current

        Args:
            None

        Returns:
            callable: slope (A/V), intercept = 0A
        """
        # return self.op_currenttovoltage().inverse
        return linop.LinearOperator(1.0 / self.Rout, units.Quantity(0, "ampere"))

    def op_fluxtocurrent(self, energy, weights=None):
        """Operator to convert flux to current

        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            callable: slope (C/ph), intercept(C/s)
        """
        Cs = self._chargepersamplephoton(energy, weights=weights, keepdims=True)
        return linop.LinearOperator(Cs, self.darkcurrent)

    def op_currenttoflux(self, energy, weights=None):
        """Operator to convert current to flux

        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            callable: slope (ph/C), intercept(ph/s)
        """
        return self.op_fluxtocurrent(energy, weights=weights).inverse

    def op_currenttocps(self):
        """Operator to convert current to counts-per-second

        Args:
            None

        Returns:
            callable: slope (cps/C), intercept(cps)
        """
        return self.oscillator.op_voltagetocps() * self.op_currenttovoltage()

    def op_cpstocurrent(self):
        """Operator to convert counts-per-second to current

        Args:
            None

        Returns:
            callable: slope (C/cps), intercept(C/s)
        """
        return self.op_voltagetocurrent() * self.oscillator.op_cpstovoltage()

    def op_countstocurrent(self, time):
        """Operator to convert counts to current

        Args:
            time: sec

        Returns:
            callable: slope (C/cts), intercept(C/s)
        """
        return self.op_cpstocurrent * self.op_countstocps(time)

    def op_voltagetocps(self):
        """Operator to convert voltage to counts-per-second

        Args:
            None

        Returns:
            callable: slope (cps/V), intercept(cps)
        """
        return self.oscillator.op_voltagetocps()

    def op_cpstovoltage(self):
        """Operator to convert counts-per-second to voltage

        Args:
            None

        Returns:
            callable: slope (V/cps), intercept(V)
        """
        return self.oscillator.op_cpstovoltage()

    def op_countstovoltage(self, time):
        """Operator to convert counts to voltage

        Args:
            None

        Returns:
            callable: slope (V/cts), intercept(V)
        """
        return self.op_cpstovoltage() * self.op_countstocps(time)

    def op_cpstoflux(self, energy, weights=None):
        """Operator to convert counts-per-second to flux

        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            callable: slope (ph/cps), intercept(ph/s)
        """
        return self.op_currenttoflux(energy, weights=weights) * self.op_cpstocurrent()

    def op_fluxtocps(self, energy, weights=None):
        """Operator to convert flux to counts-per-second

        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            callable: slope (cps/ph), intercept(cps)
        """
        return self.op_currenttocps() * self.op_fluxtocurrent(energy, weights=weights)

    def op_countstocps(self, time):
        return linop.LinearOperator(
            units.Quantity(1.0 / time, "Hz"), units.Quantity(0, "Hz")
        )

    def op_countstoflux(self, energy, time, weights=None):
        """Operator to convert counts to flux

        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            time(num): sec
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            callable: slope (ph/cts), intercept(ph/s)
        """
        return self.op_cpstoflux(energy, weights=weights) * self.op_countstocps(time)

    def op_voltagetoflux(self, energy, weights=None):
        """Operator to convert voltage to flux

        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            callable: slope (ph/s/V), intercept(ph/s)
        """
        return (
            self.op_currenttoflux(energy, weights=weights) * self.op_voltagetocurrent()
        )

    def op_fluxtovoltage(self, energy, weights=None):
        """Operator to convert flux to voltage

        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            callable: slope (V.s/ph), intercept(V)
        """
        return self.op_currenttovoltage() * self.op_fluxtocurrent(
            energy, weights=weights
        )

    def fluxtocurrent(self, energy, flux, weights=None):
        """
        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            flux(num or array): ph/s
                                shape = [nSource x] n
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            num or array: current (A)
        """
        op = self.op_fluxtocurrent(energy, weights=weights)
        return op(units.Quantity(flux, "hertz")).to("ampere")

    def currenttoflux(self, energy, current, weights=None):
        """
        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            current(num or array): A
                                   shape = [nSource x] n
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            num or array: flux (ph/s)
        """

        op = self.op_currenttoflux(energy, weights=weights)
        return op(units.Quantity(current, "ampere")).to("hertz")

    def cpstoflux(self, energy, cps, weights=None):
        """
        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            cps(num or array): cts/s
                               shape = [nSource x] n
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            num or array: flux (ph/s)
        """
        op = self.op_cpstoflux(energy, weights=weights)
        return op(units.Quantity(cps, "hertz")).to("hertz")

    def countstoflux(self, energy, time, counts, weights=None):
        """
        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            time(num): s
            cps(num or array): cts/s
                               shape = [nSource x] n
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            num or array: flux (ph/s)
        """
        op = self.op_countstoflux(energy, time, weights=weights)
        return op(units.Quantity(counts, "dimensionless")).to("hertz")

    def fluxtocps(self, energy, flux, weights=None):
        """
        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            flux(num or array): ph/s
                                shape = [nSource x] n
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            num or array: cps (cts/s)
        """
        op = self.op_fluxtocps(energy, weights=weights)
        return op(units.Quantity(flux, "hertz")).to("hertz")

    def currenttocps(self, current):
        """
        Args:
            current(num or array): A

        Returns:
            num or array: cps (cts/s)
        """
        op = self.op_currenttocps()
        return op(units.Quantity(current, "ampere")).to("hertz")

    def cpstocurrent(self, cps):
        """
        Args:
            cps(num or array): cts/s

        Returns:
            num or array: current (A)
        """
        op = self.op_cpstocurrent()
        return op(units.Quantity(cps, "hertz")).to("ampere")

    def countstocurrent(self, time, counts):
        """
        Args:
            counts(num or array): cts

        Returns:
            num or array: current (A)
        """
        op = self.op_countstocurrent(time)
        return op(units.Quantity(counts, "dimensionless")).to("ampere")

    def voltagetocps(self, voltage):
        """
        Args:
            voltage(num or array): V

        Returns:
            num or array: cps (cts/s)
        """
        op = self.op_voltagetocps()
        return op(units.Quantity(voltage, "volt")).to("hertz")

    def cpstovoltage(self, cps):
        """
        Args:
            cps(num or array): cts/s

        Returns:
            num or array: voltage (V)
        """
        op = self.op_cpstovoltage()
        return op(units.Quantity(cps, "hertz")).to("volt")

    def countstocps(self, time, counts):
        """
        Args:
            counts(num or array): cts

        Returns:
            num or array: voltage (cps)
        """
        op = self.op_countstocps(time)
        return op(units.Quantity(counts, "dimensionless")).to("Hz")

    def countstovoltage(self, time, counts):
        """
        Args:
            cps(num or array): cts

        Returns:
            num or array: voltage (V)
        """
        op = self.op_cpstovoltage(time)
        return op(units.Quantity(counts, "dimensionless")).to("volt")

    def voltagetoflux(self, energy, voltage, weights=None):
        """
        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            voltage(num or array): V
                                   shape = [nSource x] n
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            num or array: flux (ph/s)
        """
        op = self.op_voltagetoflux(energy, weights=weights)
        return op(units.Quantity(voltage, "volt")).to("hertz")

    def fluxtovoltage(self, energy, flux, weights=None):
        """
        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            flux(num or array): ph/s
                                shape = [nSource x] n
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            num or array: voltage (V)
        """
        op = self.op_fluxtovoltage(energy, weights=weights)
        return op(units.Quantity(flux, "hertz")).to("volt")

    def currenttovoltage(self, current):
        """
        Args:
            cps(num or array): A

        Returns:
            num or array: voltage (V)
        """
        op = self.op_currenttovoltage()
        return op(units.Quantity(current, "A")).to("volt")

    def fluxop(self, energy, response, weights=None, time=None):
        """Operator to convert diode response flux.

        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            response(num or array): shape = [nSource x] n
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines
            time(Optional(num)): sec

        Returns:
            op(linop): raw diode conversion operator
        """
        if time is None:
            default = "hertz"
        else:
            default = "dimensionless"

        response = units.Quantity(response, default)
        try:
            response.to("Hz")
        except pinterrors.DimensionalityError:
            try:
                response.to("dimensionless")
                if time is None:
                    raise RuntimeError("Need exposure time to convert counts to flux.")
            except pinterrors.DimensionalityError:
                try:
                    response.to("A")
                except:
                    response.to("V")
                    return self.op_voltagetoflux(energy, weights=weights)
                else:
                    return self.op_currenttoflux(energy, weights=weights)
            else:
                return self.op_countstoflux(energy, time, weights=weights)
        else:
            return self.op_cpstoflux(energy, weights=weights)

    def responsetoflux(self, energy, response, weights=None, time=None):
        """Convert diode response to flux

        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            response(num or array): cts (default when time given), cps (default when time not given), A or V
                                    shape = [nSource x] n
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines
            time(Optional(num)): sec

        Returns:
            num or array: flux (ph/s)
        """

        if time is None:
            default = "hertz"
        else:
            default = "dimensionless"

        response = units.Quantity(response, default)
        try:
            return self.cpstoflux(energy, response.to("hertz"), weights=weights)
        except pinterrors.DimensionalityError:
            try:
                tmp = response.to("dimensionless")
                if time is None:
                    raise RuntimeError("Need exposure time to convert counts to flux.")
                return self.countstoflux(energy, time, tmp, weights=weights)
            except pinterrors.DimensionalityError:
                try:
                    return self.currenttoflux(energy, response.to("A"), weights=weights)
                except:
                    return self.voltagetoflux(energy, response.to("V"), weights=weights)

    def gainfromresponse(
        self, response_after, response_before, energy=None, weights=None, time=None
    ):
        """Try to guess the diode gain

        Args:
            response_after(num): diode response after gain application (cts, cps)
            response_before(num): diode response before gain application (A, V, ph/s)
            energy(Optional(num or array)): keV
            weights(Optional(num or array)): line fractions
            time(Optional(num)): sec

        Returns:
            gain(num): V/A or A
        """

        # Voltage from the signal on which the gain is already applied
        if time is None:
            default = "Hz"
        else:
            default = "dimensionless"
        response_after = units.Quantity(response_after, default)
        try:
            V_meas = self.cpstovoltage(response_after.to("Hz"))
        except pinterrors.DimensionalityError:
            tmp = response.to("dimensionless")
            if time is None:
                raise RuntimeError("Need exposure time to convert counts to flux.")
            V_meas = self.countstovoltage(time, tmp, weights=weights)

        # Voltage from the signal before the gain is applied
        response_before = units.Quantity(response_before, "Hz")
        try:
            V_calc = self.fluxtovoltage(
                energy, response_before.to("Hz"), weights=weights
            )
        except pinterrors.DimensionalityError:
            try:
                V_calc = self.currenttovoltage(response_before.to("A"))
            except pinterrors.DimensionalityError:
                V_calc = response_before.to("V")

        if self.gain.units == ureg.Unit("A"):
            # V_meas: inverse proportional to the real gain
            # V_calc: inverse proportional to the assumed gain
            r = units.magnitude(V_calc / V_meas, "dimensionless")
        else:
            # V_meas: proportional to the real gain
            # V_calc: proportional to the assumed gain
            r = units.magnitude(V_meas / V_calc, "dimensionless")

        r = np.median(r[np.isfinite(r)])
        if np.isfinite(r and r > 0):
            return self.gainrounder(self.gain * r)
        else:
            return None

    def xrfnormop(self, energy, expotime, reference, referencetime=None, weights=None):
        """Operator to convert the raw diode signal to a flux normalizing signal.
           Usage: Inorm = I/op(iodet)

            XRF flux normalization:

            Measured XRF:
            Ixrf(cts) = F(cps).t(s).cxrf

            Desired XRF signal:
            Ixrf(cts)/norm = Fref(cps).tref(s).cxrf

            Normalization function to be apply on the raw diode signal Idiode
            norm = F.t/(Fref.tref) = cpstoflux(Idiode/t).t/(Fref.tref) = op(Idiode)
            op: x-> cpstoflux(x/t)/Fref.t/tref

            In case reference in counts instead of photons/sec
            Fref = round_sig(cpstoflux(Idioderef/t),2)

        Args:
            energy(num or array): source energies (keV)
                                  shape = [nSource x] nSourceLines
            expotime(num): sec
            reference(num): iodet (counts) or flux (photons/sec) to which the data should be normalized
            referencetime(Optional(num)): time to which the data should be normalized
            weights(Optional(num or array): source line weights
                                            shape= [nSource x] nSourceLines

        Returns:
            op(linop): raw diode conversion operator
            Fref(num): flux in photons/s to which the data is normalized after data/op(diode)
            tref(num): time in s to which the data is normalized after data/op(diode)
        """
        # Convert from counts to photons/sec
        # op: x-> cpstoflux(x/t)
        t = units.Quantity(expotime, "s")
        op = self.op_countstoflux(energy, t, weights=weights)

        # Reference flux to which the XRF signal should be normalized
        if reference.units == ureg.hertz:  # photons/sec
            Fref = reference
        elif reference.units == ureg.dimensionless:  # counts
            # Fref = op(Idioderef)
            Fref = units.Quantity(
                round_sig(units.magnitude(op(reference), "hertz"), 2), "hertz"
            )
        else:
            raise RuntimeError(
                "Reference {} should be in photons/sec (flux) or counts (iodet).".format(
                    reference
                )
            )
        # Convert from counts to counts at reference flux "Fref" and reference time "tref"
        if referencetime is not None:
            tref = units.Quantity(referencetime, "s")
            op2 = linop.LinearOperator(
                units.Quantity(t / (Fref * tref), "s"),
                units.Quantity(0, "dimensionless"),
            )
        else:
            op2 = linop.LinearOperator(
                units.Quantity(1.0 / Fref, "s"), units.Quantity(0, "dimensionless")
            )
            tref = t
        op = op2 * op
        op.m = units.magnitude(op.m, "dimensionless")
        op.b = units.magnitude(op.b, "dimensionless")
        return op, Fref.to("hertz").magnitude, tref.to("s").magnitude

    def calibratedark(self, darkresponse, time=None):
        """Dark current from response

        Args:
            darkresponse(num or array): measured dark (cts, cps, A)

        Returns:
            None
        """
        if time is None:
            default = "Hz"
        else:
            default = "dimensionless"

        def func(x):
            return np.clip(np.nanmedian(x.magnitude), 0, None)

        darkresponse = units.Quantity(darkresponse, default)
        try:
            self.darkcurrent = self.cpstocurrent(func(darkresponse.to("Hz")))
        except pinterrors.DimensionalityError:
            try:
                tmp = darkresponse.to("dimensionless")
                if time is None:
                    raise RuntimeError("Need exposure time to convert counts to flux.")
                self.darkcurrent = self.countstocurrent(time, func(tmp))
            except pinterrors.DimensionalityError:
                self.darkcurrent = func(darkresponse.to("A"))

        self.darkcurrent = np.clip(self.darkcurrent.magnitude, 0, None)

        logger.info("Calibrated dark current: {:~e}".format(self.darkcurrent.to("e/s")))

    def calibrateF0(self, darkresponse, time=None):
        """F0 from response

        Args:
            darkresponse(num or array): measured dark (cts, cps)

        Returns:
            None
        """
        if time is None:
            default = "Hz"
        else:
            default = "dimensionless"

        def func(x):
            return units.Quantity(np.clip(np.nanmedian(x.magnitude), 0, None), "Hz")

        darkresponse = units.Quantity(darkresponse, default)
        try:
            darkcps_meas = func(darkresponse.to("Hz"))
        except pinterrors.DimensionalityError:
            tmp = darkresponse.to("dimensionless")
            if time is None:
                raise RuntimeError("Need exposure time to convert counts to flux.")
            darkcps_meas = func(self.countstocps(time, tmp))

        darkcps_calc = self.currenttocps(self.darkcurrent)
        self.oscillator.F0 += darkcps_meas - darkcps_calc
        self.oscillator.F0 = np.clip(self.oscillator.F0.magnitude, 0, None)

        logger.info("Calibrated oscillator offset: {:~e}".format(self.oscillator.F0))

    def fluxcpsinfo(self, energy, weights=None):
        """
        Args:
            energy(num or array): source energies in keV (shape: nSourceLines)
            weights(num or array): source line weights (shape: nSourceLines)

        Returns:
            dict
        """
        ret = collections.OrderedDict()
        Cs = self._chargepersamplephoton(energy, weights=weights)
        ret["Energy"] = "{} keV".format(energy)
        ret["Charge/ph"] = "{:~f}".format(Cs.to("e"))
        ret["Dark"] = "{:~e}".format(self.darkcurrent.to("e/s"))
        return ret

    def propagate(self, N, energy, tframe=None, nframe=None, weights=None):
        """Error propagation of a number of photons.

        Args:
            N(unumpy.uarray): incomming number of photons with uncertainties
            energy(numpy.array): associated energies
            tframe(num|numpy.array): time per frame (sec)
            nframe(num|numpy.array): number of frames (sec)

        Returns:
            uncertainties.core.Variable or numpy.array(uncertainties.core.Variable): detector signal in DU
        """
        if tframe is None:
            ValueError("Frame exposure time not specified.")
        if nframe is None:
            ValueError("Number of frames not specified.")
        # TODO:
        Nout = None
        return Nout  # units: DU


class CalibratedPNdiode(PNdiode):
    """
    A pn-diode with a known spectral responsivity
    """

    def __init__(
        self, energy=None, response=None, model=True, fitresponse=True, **kwargs
    ):
        """
        Args:
            energy(array): keV
            response(array): spectral responsivity (A/W)
            model(Optional(bool)): use the fitted spectral responsivity or interpolate the given response
            fitresponse(Optional(bool)): modify thickness and Ehole to response
        """
        super(CalibratedPNdiode, self).__init__(**kwargs)
        self.model = model
        self.menergy = energy
        self.mresponse = response
        self.fitresponse = fitresponse
        self._init_response()

    def _init_response(self):
        energy = self.menergy
        response = self.mresponse.to("A/W")  # or e/eV
        if self.fitresponse:
            self.fit_spectral_responsivity(energy, response)
        self.finterpol = scipy.interpolate.interp1d(
            energy, response, bounds_error=False, fill_value=np.nan
        )

    def __getstate__(self):
        state = super(CalibratedPNdiode, self).__getstate__()
        state["model"] = self.model
        state["menergy"] = self.menergy
        state["mresponse"] = self.mresponse
        state["fitresponse"] = self.fitresponse
        return state

    def __setstate__(self, state):
        super(CalibratedPNdiode, self).__setstate__(state)
        self.model = state["model"]
        self.menergy = state["menergy"]
        self.mresponse = state["mresponse"]
        self.fitresponse = state["fitresponse"]
        self._init_response()

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if not (super(CalibratedPNdiode, self).__eq__(other)):
                return False
            return (
                self.model == other.model
                and all(self.menergy == other.menergy)
                and all(self.mresponse == other.mresponse)
                and self.fitresponse == other.fitresponse
            )
        else:
            return False

    def spectral_responsivity(self, energy):
        """
        Args:
            energy(num or array): keV

        Returns:
            num or array: A/W
        """
        if self.model:
            r = super(CalibratedPNdiode, self).spectral_responsivity(energy)
        else:
            r = self.finterpol(energy)
            ind = np.isnan(r)
            try:
                if any(ind):
                    r[ind] = super(CalibratedPNdiode, self).spectral_responsivity(
                        energy[ind]
                    )
            except TypeError:
                if ind:
                    r = super(CalibratedPNdiode, self).spectral_responsivity(energy)
        return units.Quantity(r, "A/W")


class NonCalibratedPNdiode(PNdiode):
    """
    A pn-diode with an unknown spectral responsivity
    """

    def __init__(self, **kwargs):
        super(NonCalibratedPNdiode, self).__init__(**kwargs)

    def calibrate(
        self,
        response,
        sampleflux,
        energy,
        weights=None,
        caliboption="optics",
        fixdark=False,
        fluxmin=0,
        fluxmax=np.inf,
    ):
        """
        Calibrate with another diode measuring the flux at the sample position

        Args:
            response(array): count rate (Hz) or current (A) measured by this diode (nResponse)
            sampleflux(array): flux measured at the sample position (Hz, nResponse)
            energy(num or array): source lines (keV, nSourceLines)
            weights(Optional(num or array)): source line weights (nSourceLines)
            caliboption(str): "optics", "solidangle" or "thickness"
            fixdark(Optional(num)): fix dark current

        Returns:
            None
        """
        # I(A) = Is(ph/s).slope + intercept
        #      = Is(ph/s).Cs + D(A)
        # Cs: charge per sample photon
        response = units.Quantity(response, "hertz")
        try:
            current = response.to("A")
        except pinterrors.DimensionalityError:
            current = self.cpstocurrent(response)

        x = units.umagnitude(sampleflux, "hertz")
        x = instance.asarray(x)
        indfit = (x >= units.umagnitude(fluxmin, "hertz")) & (
            x <= units.umagnitude(fluxmax, "hertz")
        )
        npts = sum(indfit)

        if npts < 1:
            raise RuntimeError(
                "Not enough data points with a flux in [{:e},{:e}] (data range [{:e},{:e}])".format(
                    fluxmin, fluxmax, np.min(x), np.max(x)
                )
            )
        fixdark |= npts == 1

        x = x[indfit]
        if fixdark:
            # Fit dark-subtracted current vs flux
            intercept = units.magnitude(self.darkcurrent, "ampere")
            y = units.magnitude(current, "ampere")
            y = instance.asarray(y)[indfit]
            if npts == 1:
                slope = (y[0] - intercept) / x[0]
            else:
                slope = fit1d.linfit_zerointercept(x, y - intercept)
        else:
            # Fit current vs flux
            y = units.magnitude(current, "ampere")
            y = y[indfit]
            slope, intercept = fit1d.linfit(x, y)
            if intercept < 0:
                intercept = 0
                slope = fit1d.linfit_zerointercept(x, y)

            # Set dark current:
            self.darkcurrent = intercept

        # Correlation coefficient
        ycalc = intercept + slope * x
        if npts == 1:
            R2 = np.nan
        else:
            R2 = 1 - sum((y - ycalc) ** 2) / sum((y - np.mean(y)) ** 2)

        # Set diode thickness, solid angle or transmission
        Cscalib = units.Quantity(slope, "coulomb")
        Cscalc = self._calibrate_chargepersamplephoton(
            energy, Cscalib, weights=weights, caliboption=caliboption
        )

        info = (
            "Diode calibration:\n Energy: {} keV"
            "\n Electron-hole pairs per photon hitting the sample: {:~} (experiment: {:~})"
            "\n Electron-hole pairs per second (dark): {:~e} "
            "\n R^2 = {}".format(
                energy, Cscalc.to("e"), Cscalib.to("e"), self.darkcurrent.to("e/s"), R2
            )
        )
        logger.info(info)

        ret = self.fluxcpsinfo(energy, weights=weights)
        ret["R$^2$"] = R2
        return ret


class SXM_PTB(CalibratedPNdiode):
    """PTB diode used to calibrate"""

    aliases = ["ptb"]

    def __init__(self, **kwargs):
        kwargs["attenuators"] = {}
        kwargs["attenuators"]["Detector"] = {
            "material": element.Element("Si"),
            "thickness": 30e-4,
        }
        kwargs["ehole"] = constants.eholepair_si()
        kwargs["model"] = kwargs.get("model", True)

        ptb = np.loadtxt(resource_filename("id21/ptb.dat"))
        energy = ptb[:, 0]  # keV
        response = ureg.Quantity(ptb[:, 1], "milliampere/watt")

        super(SXM_PTB, self).__init__(
            gain=ureg.Quantity(1e5, "volt/ampere"),
            gainrounder=GainRounder(base=10),
            darkcurrent=ureg.Quantity(0, "ampere"),
            energy=energy,
            response=response,
            fitresponse=True,
            beforesample=False,
            **kwargs
        )


class SXM_IDET(CalibratedPNdiode):
    """Centronic OSD 50-3T
    Keithley K428 (10V max analog output)
    NOVA N101VTF voltage-to-frequency converter (Fmax=1e6Hz, F0=0Hz, Vmax=10V)
    P201 counter board
    """

    aliases = ["idet"]

    def __init__(self, **kwargs):
        kwargs["attenuators"] = {}
        kwargs["attenuators"]["Detector"] = {
            "material": element.Element("Si"),
            "thickness": 300e-4,
        }
        kwargs["ehole"] = constants.eholepair_si()
        kwargs["model"] = kwargs.get("model", False)

        ird = np.loadtxt(resource_filename("id21/ird.dat"))

        npop = kwargs.pop("npop", None)
        if npop is None:
            npop = 4
        j = ird.shape[0] - npop
        energy = ird[:j, 0]  # keV
        responseratio = ird[:j, 1]

        energyadd = 8.4
        if energy[-1] < energyadd:
            responseratio = np.append(responseratio, 3.22)
            energy = np.append(energy, energyadd)

        absdiode = SXM_PTB(model=True)
        response = responseratio * absdiode.spectral_responsivity(energy)

        vtof = Oscillator(
            Fmax=ureg.Quantity(1e6, "Hz"),
            F0=ureg.Quantity(0, "Hz"),
            Vmax=ureg.Quantity(10, "volt"),
        )

        super(SXM_IDET, self).__init__(
            gain=ureg.Quantity(1e5, "volt/ampere"),
            gainrounder=GainRounder(base=10),
            darkcurrent=ureg.Quantity(0, "ampere"),
            energy=energy,
            response=response,
            fitresponse=False,
            beforesample=False,
            oscillator=vtof,
            **kwargs
        )


class SXM_FDET(NonCalibratedPNdiode):
    """
    Keithley K428 (10V max analog output)
    NOVA N101VTF voltage-to-frequency converter (Fmax=1e6Hz, F0=0Hz, Vmax=10V)
    P201 counter board
    """

    aliases = ["fdet"]

    def __init__(self, **kwargs):
        kwargs["attenuators"] = {}
        kwargs["attenuators"]["Detector"] = {
            "material": element.Element("Si"),
            "thickness": 0.1,
        }
        kwargs["ehole"] = constants.eholepair_si()

        vtof = Oscillator(
            Fmax=ureg.Quantity(1e6, "Hz"),
            F0=ureg.Quantity(0, "Hz"),
            Vmax=ureg.Quantity(10, "volt"),
        )

        super(SXM_FDET, self).__init__(
            gain=ureg.Quantity(1e5, "volt/ampere"),
            gainrounder=GainRounder(base=10),
            darkcurrent=ureg.Quantity(0, "ampere"),
            oscillator=vtof,
            beforesample=False,
            **kwargs
        )


class SXM_IODET1(NonCalibratedPNdiode):
    """International Radiation Detectors (IRD), AXUV-PS1-S
    Keithley K428 (10V max analog output)
    NOVA N101VTF voltage-to-frequency converter (Fmax=1e6Hz, F0=0Hz, Vmax=10V)
    P201 counter board
    """

    aliases = ["iodet1"]

    def __init__(self, **kwargs):
        kwargs2 = {}
        if "source" in kwargs:
            kwargs2["source"] = kwargs.pop("source")
        else:
            kwargs2["source"] = "synchrotron"
        if instance.isstring(kwargs2["source"]):
            kwargs2["source"] = xraysources.factory(kwargs2["source"])
        kwargs2["detector"] = self
        geometry = diodegeometries.factory("sxm_iodet1", **kwargs2)

        kwargs["attenuators"] = {}
        kwargs["attenuators"]["Detector"] = {
            "material": element.Element("Si"),
            "thickness": 0.1,
        }
        kwargs["ehole"] = constants.eholepair_si()

        window = compoundfromname.compoundfromname("silicon nitride")
        coating = element.Element("Ti")
        secondarytarget = multilayer.Multilayer(
            material=[coating, window], thickness=[500e-7, 500e-7], geometry=geometry
        )

        vtof = Oscillator(
            Fmax=ureg.Quantity(1e6, "Hz"),
            F0=ureg.Quantity(0, "Hz"),
            Vmax=ureg.Quantity(10, "volt"),
        )

        super(SXM_IODET1, self).__init__(
            gain=ureg.Quantity(1e5, "volt/ampere"),
            gainrounder=GainRounder(base=10),
            darkcurrent=ureg.Quantity(0, "ampere"),
            oscillator=vtof,
            secondarytarget=secondarytarget,
            beforesample=True,
            **kwargs
        )


class SXM_IODET2(NonCalibratedPNdiode):
    """International Radiation Detectors (IRD), AXUV-PS1-S
    Keithley K428 (10V max analog output)
    NOVA N101VTF voltage-to-frequency converter (Fmax=1e6Hz, Vmax=0Hz, Vmax=10V)
    P201 counter board
    """

    aliases = ["iodet2"]

    def __init__(self, **kwargs):
        kwargs2 = {}
        if "source" in kwargs:
            kwargs2["source"] = kwargs.pop("source")
        else:
            kwargs2["source"] = "synchrotron"
        if instance.isstring(kwargs2["source"]):
            kwargs2["source"] = xraysources.factory(kwargs2["source"])
        kwargs2["detector"] = self
        geometry = diodegeometries.factory("sxm_iodet2", **kwargs2)

        kwargs["attenuators"] = {}
        kwargs["attenuators"]["Detector"] = {
            "material": element.Element("Si"),
            "thickness": 0.1,
        }
        kwargs["ehole"] = constants.eholepair_si()

        window = compoundfromname.compoundfromname("silicon nitride")
        secondarytarget = multilayer.Multilayer(
            material=[window], thickness=[500e-7], geometry=geometry
        )

        vtof = Oscillator(
            Fmax=ureg.Quantity(1e6, "Hz"),
            F0=ureg.Quantity(0, "Hz"),
            Vmax=ureg.Quantity(10, "volt"),
        )

        super(SXM_IODET2, self).__init__(
            gain=ureg.Quantity(1e5, "volt/ampere"),
            gainrounder=GainRounder(base=10),
            darkcurrent=ureg.Quantity(0, "ampere"),
            oscillator=vtof,
            secondarytarget=secondarytarget,
            beforesample=True,
            **kwargs
        )


class XRD_IDET(NonCalibratedPNdiode):
    aliases = ["microdiff_pico1"]

    def __init__(self, **kwargs):
        kwargs["attenuators"] = {}
        kwargs["attenuators"]["Detector"] = {
            "material": element.Element("Si"),
            "thickness": 0.1,
        }
        kwargs["ehole"] = constants.eholepair_si()
        super(XRD_IDET, self).__init__(
            gain=ureg.Quantity(2.1e-6, "ampere"),
            gainrounder=GainRounder(base=10, m=2.1),
            darkcurrent=ureg.Quantity(0, "ampere"),
            beforesample=False,
            **kwargs
        )


class ID16B_IT(NonCalibratedPNdiode):
    """
    Keithley K428 (2V max analog output)
    V2F100 voltage-to-frequency converter (Fmax=50e6Hz, F0=?Hz, Vmax=2.5V)
    P201 counter board
    """

    aliases = ["id16b_It"]

    def __init__(self, **kwargs):
        kwargs["attenuators"] = {}
        kwargs["attenuators"][self.ATMOSPHERELABEL] = {
            "material": compoundfromname.compoundfromname("air"),
            "thickness": 80.0,
        }
        kwargs["attenuators"]["Detector"] = {
            "material": element.Element("Si"),
            "thickness": 500e-4,
        }
        kwargs["ehole"] = constants.eholepair_si()

        vtof = Oscillator(
            Fmax=ureg.Quantity(50e6, "Hz"),
            F0=ureg.Quantity(0, "Hz"),
            Vmax=ureg.Quantity(2.5, "volt"),
        )

        super(ID16B_IT, self).__init__(
            gain=ureg.Quantity(2.1e-6, "ampere"),
            gainrounder=GainRounder(base=10, m=2.1),
            darkcurrent=ureg.Quantity(0, "ampere"),
            oscillator=vtof,
            Vmax=ureg.Quantity(2.1, "volt"),
            beforesample=False,
            **kwargs
        )


class ID16B_I0(NonCalibratedPNdiode):
    """
    Keithley K428 (2V max analog output)
    V2F100 voltage-to-frequency converter (Fmax=50e6Hz, F0=?Hz, Vmax=2.5V)
    P201 counter board
    """

    aliases = ["id16b_I0"]

    def __init__(self, **kwargs):

        kwargs["attenuators"] = {}
        kwargs["attenuators"]["Detector"] = {
            "material": element.Element("Si"),
            "thickness": 500e-4,
        }
        kwargs["ehole"] = constants.eholepair_si()

        vtof = Oscillator(
            Fmax=ureg.Quantity(50e6, "Hz"),
            F0=ureg.Quantity(0, "Hz"),
            Vmax=ureg.Quantity(2.5, "volt"),
        )

        super(ID16B_I0, self).__init__(
            gain=ureg.Quantity(2.1e-6, "ampere"),
            gainrounder=GainRounder(base=10, m=2.1),
            darkcurrent=ureg.Quantity(0, "ampere"),
            oscillator=vtof,
            Vmax=ureg.Quantity(2.1, "volt"),
            beforesample=True,
            **kwargs
        )


class ID16B_IC(NonCalibratedPNdiode):
    """
    Keithley K428 (2V max analog output)
    V2F100 voltage-to-frequency converter (Fmax=50e6Hz, F0=?Hz, Vmax=2.5V)
    P201 counter board
    """

    aliases = ["id16b_IC"]

    def __init__(self, **kwargs):
        kwargs2 = {}
        if "source" in kwargs:
            kwargs2["source"] = kwargs.pop("source")
        else:
            kwargs2["source"] = "synchrotron"
        if instance.isstring(kwargs2["source"]):
            kwargs2["source"] = xraysources.factory(kwargs2["source"])
        kwargs2["detector"] = self
        geometry = diodegeometries.factory("id16b_ic", **kwargs2)

        kwargs["attenuators"] = {}
        kwargs["attenuators"]["Detector"] = {
            "material": element.Element("Si"),
            "thickness": 500e-4,
        }
        kwargs["ehole"] = constants.eholepair_si()

        window = compoundfromname.compoundfromname(
            "silicon nitride"
        )  # TODO: mirror material + thickness???
        secondarytarget = multilayer.Multilayer(
            material=[window], thickness=[500e-7], geometry=geometry
        )

        vtof = Oscillator(
            Fmax=ureg.Quantity(50e6, "Hz"),
            F0=ureg.Quantity(0, "Hz"),
            Vmax=ureg.Quantity(2.5, "volt"),
        )

        super(ID16B_IC, self).__init__(
            gain=ureg.Quantity(2.1e-6, "ampere"),
            gainrounder=GainRounder(base=10, m=2.1),
            darkcurrent=ureg.Quantity(0, "ampere"),
            oscillator=vtof,
            Vmax=ureg.Quantity(2.1, "volt"),
            secondarytarget=secondarytarget,
            beforesample=True,
            **kwargs
        )


factory = PNdiode.factory
registry = PNdiode.clsregistry

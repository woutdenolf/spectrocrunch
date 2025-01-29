# -*- coding: utf-8 -*-

from . import base
from ..materials import compoundfromname
from ..materials import element
from ..utils import constants
from ..utils.classfactory import with_metaclass
from ..utils import instance
from ..utils import units
from ..math import linalg

import scipy.special
import numpy as np
from copy import copy


class XRFDetector(with_metaclass(base.CentricCone)):
    """
    Line profile:

        1. R. van Grieken, A. Markowicz: Handbook of X-ray Spectrometry, CRC Press, Boca Raton (2001)
           peak, tail and step have fixed area ratio's
           area includes peak, tail and step
        2. PyMca
           peak and tail have a fixed area ratio
           peak and step have a fixed height ratio
           area includes only peak
    """

    def __init__(
        self,
        mcazero=None,
        mcagain=None,
        mcanoise=None,
        mcafano=None,
        shape_fixedarearatios=None,
        shape_pymca=None,
        shape_conversionenergy=None,
        **kwargs
    ):
        """
        Args:
            mcazero(num):
            mcagain(num):
            mcanoise(num):
            mcafano(num):
            shape_fixedarearatios(Optional(dict)): see van Grieken
            shape_pymca(Optional(dict)): see PyMca
        """
        super(XRFDetector, self).__init__(**kwargs)

        self.mcazero = mcazero  # keV
        self.mcagain = mcagain  # keV/bin
        self.mcanoise = mcanoise  # FWHM in keV
        self.mcafano = mcafano

        self.shape_conversionenergy = shape_conversionenergy
        self.fixedarearatios = shape_fixedarearatios is not None
        if self.fixedarearatios:
            self.bpeak = shape_fixedarearatios["bpeak"]
            self.bstail = shape_fixedarearatios["bstail"]
            self.bltail = shape_fixedarearatios["bltail"]
            self.bstep = shape_fixedarearatios["bstep"]
            self.stailbroadening = shape_fixedarearatios["stailbroadening"]
            self.ltailbroadening = shape_fixedarearatios["ltailbroadening"]
            self.fractions = (
                shape_fixedarearatios["stailfraction"],
                shape_fixedarearatios["ltailfraction"],
                shape_fixedarearatios["stepfraction"],
            )
        else:
            if shape_pymca is None:
                shape_pymca = {}
            self.bpeak = shape_pymca.get("bpeak", True)
            self.bstail = shape_pymca.get("bstail", False)
            self.bltail = shape_pymca.get("bltail", False)
            self.bstep = shape_pymca.get("bstep", False)
            self.stailslope_ratio = shape_pymca.get("stailslope_ratio", 0.5)
            self.ltailslope_ratio = shape_pymca.get("ltailslope_ratio", 10.0)
            self.ratios = (
                shape_pymca.get("stailarea_ratio", 0.05),
                shape_pymca.get("ltailarea_ratio", 0.02),
                shape_pymca.get("stepheight_ratio", 0.0001),
            )

    def convert(self, shape_conversionenergy=None, inplace=False):
        """Convert from a normalized peak area (including tails and step) to
        a normalized Gaussian peak area (excluding tails and step). Note that
        this conversion is only valid for one energy.
        """
        if inplace:
            o = self
        else:
            o = copy(self)
        if shape_conversionenergy is None:
            shape_conversionenergy = o.shape_conversionenergy
        else:
            o.shape_conversionenergy = shape_conversionenergy
        if shape_conversionenergy is None:
            raise RuntimeError("Specify conversion energy")
        if o.fixedarearatios:
            o.stailbroadening = o.stailbroadening
            o.ltailbroadening = o.ltailbroadening
            o.fractions = o.fractions
        else:
            o.stailslope_ratio = o.stailslope_ratio
            o.ltailslope_ratio = o.ltailslope_ratio
            o.ratios = o.ratios
        o.fixedarearatios = not o.fixedarearatios
        return o

    def __getstate__(self):
        state = super(XRFDetector, self).__getstate__()
        state["mcazero"] = self.mcazero
        state["mcagain"] = self.mcagain
        state["mcanoise"] = self.mcanoise
        state["mcafano"] = self.mcafano
        state["shape_conversionenergy"] = self.shape_conversionenergy
        state["fixedarearatios"] = self.fixedarearatios
        state["bpeak"] = self.bpeak
        state["bstail"] = self.bstail
        state["bltail"] = self.bltail
        state["bstep"] = self.bstep
        if self.fixedarearatios:
            state["stailbroadening"] = self.stailbroadening
            state["ltailbroadening"] = self.ltailbroadening
            state["fractions"] = self.fractions
        else:
            state["stailslope_ratio"] = self.stailslope_ratio
            state["ltailslope_ratio"] = self.ltailslope_ratio
            state["ratios"] = self.ratios
        return state

    def __setstate__(self, state):
        super(XRFDetector, self).__setstate__(state)
        self.mcazero = state["mcazero"]
        self.mcagain = state["mcagain"]
        self.mcanoise = state["mcanoise"]
        self.mcafano = state["mcafano"]
        self.shape_conversionenergy = state["shape_conversionenergy"]
        self.fixedarearatios = state["fixedarearatios"]
        self.bpeak = state["bpeak"]
        self.bstail = state["bstail"]
        self.bltail = state["bltail"]
        self.bstep = state["bstep"]
        if self.fixedarearatios:
            self.stailbroadening = state["stailbroadening"]
            self.ltailbroadening = state["ltailbroadening"]
            self.fractions = state["fractions"]
        else:
            self.stailslope_ratio = state["stailslope_ratio"]
            self.ltailslope_ratio = state["ltailslope_ratio"]
            self.ratios = state["ratios"]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if not (super(XRFDetector, self).__eq__(other)):
                return False
            if self.fixedarearatios != other.fixedarearatios:
                return False
            if self.fixedarearatios:
                if (
                    self.stailbroadening != other.stailbroadening
                    or self.ltailbroadening != other.ltailbroadening
                    or self.fractions != other.fractions
                ):
                    return False
            else:
                if (
                    self.stailslope_ratio != other.stailslope_ratio
                    or self.ltailslope_ratio != other.ltailslope_ratio
                    or self.ratios != other.ratios
                ):
                    return False
            return (
                self.mcazero == other.mcazero
                and self.mcagain == other.mcagain
                and self.mcanoise == other.mcanoise
                and self.mcafano == other.mcafano
                and self.shape_conversionenergy == other.shape_conversionenergy
                and self.bpeak == other.bpeak
                and self.bstail == other.bstail
                and self.bltail == other.bltail
                and self.bstep == other.bstep
            )
        else:
            return False

    @property
    def _calc_fixedarearatios(self):
        return self.fixedarearatios and self.shape_conversionenergy is not None

    @property
    def _calc_notfixedarearatios(self):
        return not self.fixedarearatios and self.shape_conversionenergy is not None

    def _calc_energy(self, energy):
        if energy is None:
            energy = self.shape_conversionenergy
        return energy

    @property
    def stailbroadening(self):
        return self._stailbroadening

    @stailbroadening.setter
    def stailbroadening(self, value):
        self._stailbroadening = value
        if self._calc_fixedarearatios:
            self.stailslope_ratio = self._calc_tail_slope_ratio(self.stailbroadening)

    @property
    def ltailbroadening(self):
        return self._ltailbroadening

    @ltailbroadening.setter
    def ltailbroadening(self, value):
        self._ltailbroadening = value
        if self._calc_fixedarearatios:
            self.ltailslope_ratio = self._calc_tail_slope_ratio(self.ltailbroadening)

    @property
    def stailslope_ratio(self):
        return self._stailslope_ratio

    @stailslope_ratio.setter
    def stailslope_ratio(self, value):
        self._stailslope_ratio = value
        if self._calc_notfixedarearatios:
            self.stailbroadening = self._calc_tail_boardening(self.stailslope_ratio)

    @property
    def ltailslope_ratio(self):
        return self._ltailslope_ratio

    @ltailslope_ratio.setter
    def ltailslope_ratio(self, value):
        self._ltailslope_ratio = value
        if self._calc_notfixedarearatios:
            self.ltailbroadening = self._calc_tail_boardening(self.ltailslope_ratio)

    def _calc_tail_slope_ratio(self, broadening, energy=None):
        energy = self._calc_energy(energy)
        gsigma = np.sqrt(self.gaussianVAR(energy))
        return self._fcalc_tail_slope_ratio(broadening, gsigma)

    def _calc_tail_boardening(self, slope_ratio, energy=None):
        energy = self._calc_energy(energy)
        gsigma = np.sqrt(self.gaussianVAR(energy))
        return self._fcalc_tail_boardening(slope_ratio, gsigma)

    def _fcalc_tail_boardening(self, slope_ratio, gsigma):
        return slope_ratio / gsigma

    def _fcalc_tail_slope_ratio(self, broadening, gsigma):
        return broadening * gsigma

    @staticmethod
    def wpeak(wstail, wltail, wstep):
        s = wstail + wltail + wstep
        if s < -0.0001 or s > 1.0001:
            raise RuntimeError(
                "Sum of step and tail fractions must be between 0 and 1 ({})".format(s)
            )
        return min(max(1 - s, 0), 1)

    @property
    def fractions(self):
        wstail = self._stailfraction
        wltail = self._ltailfraction
        wstep = self._stepfraction
        return wstail, wltail, wstep

    @property
    def all_fractions(self):
        wstail = self._stailfraction
        wltail = self._ltailfraction
        wstep = self._stepfraction
        wpeak = self.wpeak(wstail, wltail, wstep)
        return wpeak, wstail, wltail, wstep

    @fractions.setter
    def fractions(self, value):
        wstail, wltail, wstep = value
        self._stailfraction = wstail
        self._ltailfraction = wltail
        self._stepfraction = wstep
        if self._calc_fixedarearatios:
            self.ratios = self._calc_ratios()

    @property
    def ratios(self):
        return self._stailarea_ratio, self._ltailarea_ratio, self._stepheight_ratio

    @ratios.setter
    def ratios(self, value):
        self._stailarea_ratio, self._ltailarea_ratio, self._stepheight_ratio = value
        if self._calc_notfixedarearatios:
            self.fractions = self._calc_fractions()

    def _calc_ratios(self, energy=None):
        # if wpeak>0:
        #   rstail = wstail/wpeak*ctail
        #   rltail = wltail/wpeak*ctail
        #   rstep = wstep/wpeak*cstep
        # else:
        #   rstail = wstail*ctail
        #   rltail = wltail*ctail
        #   rstep = wstep*cstep

        u = self._calc_energy(energy)
        gvar = self.gaussianVAR(u)
        a = 2 * gvar
        b = np.sqrt(a)
        gnorm = self._gnorm(a)
        snorm = self._snorm(u, a, b)
        gsigma = np.sqrt(gvar)
        _str = self._fcalc_tail_slope_ratio(self.stailbroadening, gsigma)
        ltr = self._fcalc_tail_slope_ratio(self.ltailbroadening, gsigma)
        stnorm = self._tnorm(u, a, b, _str)
        ltnorm = self._tnorm(u, a, b, ltr)

        return self._fcalc_ratios_single(
            self.all_fractions, _str / stnorm, ltr / ltnorm, gnorm / snorm
        )

    def _calc_fractions(self, energy=None):
        # if wpeak>0:
        #   rstail = wstail/wpeak*cstail
        #   rltail = wltail/wpeak*cltail
        #   rstep = wstep/wpeak*cstep
        # else:
        #   rstail = wstail*cstail
        #   rltail = wltail*cltail
        #   rstep = wstep*cstep

        u = self._calc_energy(energy)
        gvar = self.gaussianVAR(u)
        a = 2 * gvar
        b = np.sqrt(a)
        gnorm = self._gnorm(a)
        snorm = self._snorm(u, a, b)
        _str = self.stailslope_ratio
        ltr = self.ltailslope_ratio
        stnorm = self._tnorm(u, a, b, _str)
        ltnorm = self._tnorm(u, a, b, ltr)

        return self._fcalc_fractions_single(
            self.ratios, _str / stnorm, ltr / ltnorm, gnorm / snorm
        )

    def _fcalc_ratios_single(self, fractions, cstail, cltail, cstep):
        wpeak, wstail, wltail, wstep = fractions
        if wpeak == 0:
            wpeak = 1
        rstail = wstail / wpeak * cstail
        rltail = wltail / wpeak * cltail
        rstep = wstep / wpeak * cstep
        return rstail, rltail, rstep

    def _fcalc_ratios(self, fractions, cstail, cltail, cstep):
        if instance.isarray(cstail):
            rstail, rltail, rstep = zip(
                *tuple(
                    self._fcalc_ratios_single(fractions, cstaili, cltaili, cstepi)
                    for cstaili, cltaili, cstepi in zip(
                        cstail.flat, cltail.flat, cstep.flat
                    )
                )
            )
            rstail = np.asarray(rstail).reshape(cstail.shape)
            rltail = np.asarray(rltail).reshape(cltail.shape)
            rstep = np.asarray(rstep).reshape(cstep.shape)
        else:
            rstail, rltail, rstep = self._fcalc_ratios_single(
                fractions, cstail, cltail, cstep
            )
        return rstail, rltail, rstep

    def _fcalc_fractions_single(self, ratios, cstail, cltail, cstep):
        #   rstail = (rstail+cstail)*wstail + rstail         *wltail + rstail       *wstep
        #   rltail = rltail         *wstail + (rltail+cltail)*wltail + rltail       *wstep
        #   rstep  = rstep          *wstail + rstep          *wltail + (rstep*cstep)*wstep

        rstail, rltail, rstep = ratios
        if self.bpeak:
            A = np.array(
                [
                    [rstail + cstail, rstail, rstail],
                    [rltail, rltail + cltail, rltail],
                    [rstep, rstep, rstep + cstep],
                ]
            )
            b = np.array([rstail, rltail, rstep])
            wstail, wltail, wstep = linalg.cramer(A, b)
        else:
            wstail = rstail / cstail
            wltail = rltail / cltail
            wstep = rstep / cstep
        return wstail, wltail, wstep

    def _fcalc_fractions(self, ratios, cstail, cltail, cstep):
        if instance.isarray(cstail):
            wstail, wltail, wstep = zip(
                *tuple(
                    self._fcalc_fractions_single(ratios, cstaili, cltaili, cstepi)
                    for cstaili, cltaili, cstepi in zip(
                        cstail.flat, cltail.flat, cstep.flat
                    )
                )
            )
            wstail = np.asarray(wstail).reshape(cstail.shape)
            wltail = np.asarray(wltail).reshape(cltail.shape)
            wstep = np.asarray(wstep).reshape(cstep.shape)
        else:
            wstail, wltail, wstep = self._fcalc_fractions_single(
                ratios, cstail, cltail, cstep
            )
        return wstail, wltail, wstep

    def __str__(self):
        if self.fixedarearatios:
            wpeak, wstail, wltail, wstep = self.all_fractions
        else:
            rstail, rltail, rstep = self.ratios

        shape = ""

        if self.fixedarearatios:
            shape = "{}\n Peak fraction = {}".format(shape, wpeak)

        if self.bstail:
            if self.fixedarearatios:
                shape = (
                    "{}\n Short-tail broadening = {}"
                    "\n Short-tail fraction = {}".format(
                        shape, self.stailbroadening, wstail
                    )
                )
            else:
                shape = (
                    "{}\n Short-tail slope ratio = {}"
                    "\n Short-tail area ratio = {}".format(
                        shape, self.stailslope_ratio, rstail
                    )
                )

        if self.bltail:
            if self.fixedarearatios:
                shape = (
                    "{}\n Long-tail broadening = {}"
                    "\n Long-tail fraction = {}".format(
                        shape, self.ltailbroadening, wltail
                    )
                )
            else:
                shape = (
                    "{}\n Long-tail slope ratio = {}"
                    "\n Long-tail area ratio = {}".format(
                        shape, self.ltailslope_ratio, rltail
                    )
                )

        if self.bstep:
            if self.fixedarearatios:
                shape = "{}\n Step fraction = {}".format(shape, wstep)
            else:
                shape = "{}\n Step height ratio = {}".format(shape, rstep)

        return (
            "XRF detector:\n{}\n"
            "MCA:\n zero = {} eV\n"
            " gain = {} eV\n"
            " noise = {} eV (FWHM)\n"
            " fano = {}{}".format(
                super(XRFDetector, self).__str__(),
                self.mcazero * 1000,
                self.mcagain * 1000,
                self.mcanoise * 1000,
                self.mcafano,
                shape,
            )
        )

    def FWHMtoVAR(self, FWHM):
        return FWHM**2 / (8 * np.log(2))

    def VARtoFWHM(self, var):
        return np.sqrt(var * (8 * np.log(2)))

    def gaussianVAR(self, energy):
        """Gaussian variance (keV^2)"""
        return (
            self.FWHMtoVAR(self.mcanoise)
            + self.mcafano * self.ehole.to("keV").magnitude * energy
        )

    def gaussianFWHM(self, energy):
        """Gaussian FWHM (keV)"""
        return self.VARtoFWHM(self.gaussianVAR(energy))

    def voigtFWHM(self, energy, linewidth=0):
        fG = self.gaussianFWHM(energy)  # detection FWHM
        fL = linewidth  # transition FWHM
        return 0.5346 * fL + np.sqrt(0.2166 * fL**2 + fG**2)

    def voigtVAR(self, energy, linewidth=0):
        return self.FWHMtoVAR(self.voigtFWHM(energy, linewidth=linewidth))

    def _gnorm(self, a):
        # Normalize H.Gaussian in -inf<x<inf
        return np.sqrt(np.pi * a)

    def _snorm(self, u, a, b, approx=False):
        # Normalize H.Step in 0<=x<inf
        if approx:
            snorm = float(u)
        else:
            zero = (
                a * np.exp(-(u**2) / a) / (2 * np.sqrt(np.pi))
                - u * scipy.special.erfc(u / b) / 2.0
            )
            snorm = u + zero
        return snorm

    def _tnorm(self, u, a, b, tr, approx=False):
        # Normalize H.Tail in -inf<x<inf
        if approx:
            return float(tr)
        else:
            minusone = np.exp(a / (4.0 * tr**2) - u / tr) * scipy.special.erfc(
                (a / (2.0 * tr) - u) / b
            ) + scipy.special.erf(-u / b)
            return tr * (1 - minusone) / 2.0

    def lineprofile(
        self, x, u, linewidth=0, normalized=None, decomposed=False, onlyheight=False
    ):
        """
        Args:
            x(num|array): energies (keV, nx)
            u(num|array): peak energies (keV, nu)
            linewidth(Optional(num|array)): natural line widths
            onlyheight(Optional(bool)): line heights instead of full profile

        Returns:
            array: nx x nu
        """
        x = instance.asarray(x)[:, np.newaxis]
        u = instance.asarray(u)[np.newaxis, :]

        gvar = self.gaussianVAR(u)
        a = 2 * gvar  # 2.sigma^2
        b = np.sqrt(a)  # sqrt(2).sigma
        gnorm = self._gnorm(a)

        stail_H = 0
        ltail_H = 0
        stailslope_ratio = 0
        ltailslope_ratio = 0
        if self.fixedarearatios:
            if normalized is None:
                normalized = True

            wpeak, wstail, wltail, wstep = self.all_fractions
            stailbroadening = self.stailbroadening
            ltailbroadening = self.ltailbroadening

            bpeak = wpeak > 0 and self.bpeak
            bstail = stailbroadening > 0 and wstail > 0 and self.bstail
            bltail = ltailbroadening > 0 and wltail > 0 and self.bltail
            bstep = wstep > 0 and self.bstep

            if normalized:
                if bpeak:
                    peak_H = wpeak / gnorm
                if bstep:
                    snorm = self._snorm(u, a, b)
                    step_H = wstep / snorm
                if bstail or bltail:
                    gsigma = np.sqrt(gvar)

                    if bstail:
                        stailslope_ratio = self._fcalc_tail_slope_ratio(
                            stailbroadening, gsigma
                        )
                        stnorm = self._tnorm(u, a, b, stailslope_ratio)
                        stail_H = wstail / stnorm

                    if bltail:
                        ltailslope_ratio = self._fcalc_tail_slope_ratio(
                            ltailbroadening, gsigma
                        )
                        ltnorm = self._tnorm(u, a, b, ltailslope_ratio)
                        ltail_H = wltail / ltnorm
            else:
                bpeak = self.bpeak
                if bpeak:
                    peak_H = 1 / gnorm

                if bstep or bstail or bltail:
                    gsigma = np.sqrt(gvar)
                    stailslope_ratio = self._fcalc_tail_slope_ratio(
                        stailbroadening, gsigma
                    )
                    ltailslope_ratio = self._fcalc_tail_slope_ratio(
                        ltailbroadening, gsigma
                    )

                    snorm = self._snorm(u, a, b)
                    stnorm = self._tnorm(u, a, b, stailslope_ratio)
                    ltnorm = self._tnorm(u, a, b, ltailslope_ratio)

                    (
                        stailarea_ratio,
                        ltailarea_ratio,
                        stepheight_ratio,
                    ) = self._fcalc_ratios(
                        (wpeak, wstail, wltail, wstep),
                        stailslope_ratio / stnorm,
                        ltailslope_ratio / ltnorm,
                        gnorm / snorm,
                    )

                    if bstail:
                        stail_H = stailarea_ratio / stailslope_ratio
                    if bltail:
                        ltail_H = ltailarea_ratio / ltailslope_ratio
                    if bstep:
                        step_H = stepheight_ratio / gnorm

        else:
            if normalized is None:
                normalized = False

            stailarea_ratio, ltailarea_ratio, stepheight_ratio = self.ratios
            stailslope_ratio = self.stailslope_ratio
            ltailslope_ratio = self.ltailslope_ratio

            bpeak = self.bpeak
            bstail = stailslope_ratio > 0 and stailarea_ratio > 0 and self.bstail
            bltail = ltailslope_ratio > 0 and ltailarea_ratio > 0 and self.bltail
            bstep = stepheight_ratio > 0 and self.bstep

            if normalized:
                snorm = self._snorm(u, a, b)
                stnorm = self._tnorm(u, a, b, stailslope_ratio)
                ltnorm = self._tnorm(u, a, b, ltailslope_ratio)

                wstail, wltail, wstep = self._fcalc_fractions(
                    (stailarea_ratio, ltailarea_ratio, stepheight_ratio),
                    stailslope_ratio / stnorm,
                    ltailslope_ratio / ltnorm,
                    gnorm / snorm,
                )
                wpeak = self.wpeak(wstail, wltail, wstep)

                if bpeak:
                    peak_H = wpeak / gnorm
                if bstail:
                    stail_H = wstail / stnorm
                if bltail:
                    ltail_H = wltail / ltnorm
                if bstep:
                    step_H = wstep / snorm
            else:
                bpeak = self.bpeak
                peak_H = 1 / gnorm
                if bstail:
                    stail_H = stailarea_ratio / stailslope_ratio
                if bltail:
                    ltail_H = ltailarea_ratio / ltailslope_ratio
                if bstep:
                    step_H = stepheight_ratio / gnorm

        if not onlyheight:
            diff = x - u
            if bstep or bstail or bltail:
                argstep = diff / b  # (x-u)/(sqrt(2).sigma)

        # XRFDetector response:
        #   Gaussian: H*exp(-(x-u)^2/(2.gvar))
        #   Lorentz:  2/(pi.W)/(1+4/W^2.(x-u)^2)    (W:FWHM)
        #             W/(2.pi)/(W^2/4+(x-u)^2)
        #             y/pi/(y^2+(x-u)^2)            (y:W/2)
        #   Voigt:    H.Re(w[(x-u+i.y)/(2.sqrt(gvar))])
        #
        #   VanGrieken: H = wpeak/gnorm            (normalize in ]-inf,inf[)
        #               gnorm = sqrt(2.pi.gvar)
        #   Pymca: H = garea/gnorm
        #               gnorm = sqrt(2.pi.gvar)
        #
        #   => garea = wpeak
        if onlyheight:
            if bpeak:
                yg = peak_H
            else:
                yg = np.zeros_like(u)
        else:
            if bpeak:
                W = instance.asarray(linewidth)[np.newaxis, :]
                if W.any():
                    yg = peak_H * np.real(scipy.special.wofz((diff + 0.5j * W) / b))
                else:
                    yg = peak_H * np.exp(-(diff**2) / a)
            else:
                yg = np.zeros_like(diff)

        # Incomplete charge collection (step):
        #   Gaussian step: H.step
        #                   step = erfc[(x-u)/sqrt(2.gvar)]/2
        #
        #   VanGrieken: H = wstep/snorm                   (normalize in [0,inf[, snorm≃u)
        #   Pymca: H = stepheight_ratio*garea/gnorm
        #
        #   => stepheight_ratio = wstep/wpeak*gnorm/snorm
        if onlyheight:
            if bstep:
                ys = step_H / 2.0
            else:
                ys = np.zeros_like(u)
        else:
            if bstep:
                ys = step_H / 2.0 * scipy.special.erfc(argstep)
            else:
                ys = np.zeros_like(diff)

        # Incomplete charge collection (tail):
        #   Gaussian tail: H.tail
        #                  tail = exp[(x-u)/tr+gvar/(2.tr**2)].erfc[(x-u)/sqrt(2.gvar)+sqrt(gvar)/(sqrt(2).tr)]/2
        #                  tr = tb.sqrt(gvar)
        #
        #   VanGrieken: H = wtail/tnorm               (normalize in ]-inf,inf[, tnorm≃tr)
        #   Pymca: H = garea*tailarea_ratio/tr
        #
        #   => tailarea_ratio = wtail/wpeak*tr/tnorm
        yt = []
        for btail, tail_H, tailslope_ratio in zip(
            [bstail, bltail], [stail_H, ltail_H], [stailslope_ratio, ltailslope_ratio]
        ):
            if onlyheight:
                if btail:
                    with np.errstate(over="ignore"):
                        mexp = tail_H / 2.0 * np.exp(gvar / (2.0 * tailslope_ratio**2))
                        ind = np.isinf(mexp)
                        if ind.any():
                            mexp[ind] = 0  # lim_x->inf exp(x)*erfc(x) = 0

                    yti = mexp * scipy.special.erfc(b / (2.0 * tailslope_ratio))
                else:
                    yti = np.zeros_like(u)
            else:
                if btail:
                    with np.errstate(over="ignore"):
                        mexp = (
                            tail_H
                            / 2.0
                            * np.exp(
                                diff / tailslope_ratio
                                + gvar / (2.0 * tailslope_ratio**2)
                            )
                        )
                        ind = np.isinf(mexp)
                        if ind.any():
                            mexp[ind] = 0  # lim_x->inf exp(x)*erfc(x) = 0

                    yti = mexp * scipy.special.erfc(
                        argstep + b / (2.0 * tailslope_ratio)
                    )
                else:
                    yti = np.zeros_like(diff)
            yt.append(yti)
        yst, ylt = yt

        if decomposed:
            return yg, yst, ylt, ys
        else:
            return yg + yst + ylt + ys

    def addtopymca(self, setup, cfg):
        super(XRFDetector, self).addtopymca(setup, cfg)

        mcazero = self.mcazero
        mcagain = self.mcagain
        mcanoise = self.mcanoise
        # Compensate fano for difference in Ehole with pymca's 3.85 eV (ClassMcaTheory)
        mcafano = (
            self.mcafano
            * (self.ehole / units.Quantity(3.85, "eV")).to("dimensionless").magnitude
        )
        cfg["detector"]["zero"] = mcazero
        cfg["detector"]["gain"] = mcagain
        cfg["detector"]["noise"] = mcanoise
        cfg["detector"]["fano"] = mcafano

        # No one-to-one correspondance
        self.shape_conversionenergy = np.max(setup.energy)

        rstail, rltail, rstep = self.ratios
        cfg["peakshape"]["st_arearatio"] = rstail
        cfg["peakshape"]["st_sloperatio"] = self.stailslope_ratio
        cfg["peakshape"]["lt_arearatio"] = rltail
        cfg["peakshape"]["lt_sloperatio"] = self.ltailslope_ratio
        cfg["peakshape"]["step_heightratio"] = rstep

        bhypermet = True
        bstail = self.bstail
        bltail = self.bltail
        bstep = self.bstep
        cfg["fit"]["hypermetflag"] = 1 * bhypermet | 2 * bstail | 4 * bltail | 8 * bstep

        xmin, xmax = self.channellimits(setup)
        cfg["fit"]["xmin"] = xmin
        cfg["fit"]["xmax"] = xmax
        cfg["fit"]["use_limit"] = 1

    def channellimits(self, setup):
        energy, func = instance.asarrayf([setup.emin, setup.emax])
        return func(
            np.clip(
                np.round((energy - self.mcazero) / self.mcagain).astype(int), 0, None
            )
        )

    def loadfrompymca(self, setup, cfg):
        super(XRFDetector, self).loadfrompymca(setup, cfg)

        self.mcazero = cfg["detector"]["zero"]
        self.mcagain = cfg["detector"]["gain"]
        self.mcanoise = cfg["detector"]["noise"]
        # Compensate fano for difference in Ehole with pymca's 3.85 eV (ClassMcaTheory)
        self.mcafano = (
            cfg["detector"]["fano"]
            * (units.Quantity(3.85, "eV") / self.ehole).to("dimensionless").magnitude
        )

        # No one-to-one correspondance
        self.shape_conversionenergy = np.max(setup.energy)

        b = cfg["fit"]["hypermetflag"]
        self.bstail = (b & 2) == 2
        self.bltail = (b & 4) == 4
        self.bstep = (b & 8) == 8

        rstail = cfg["peakshape"]["st_arearatio"]
        self.stailslope_ratio = cfg["peakshape"]["st_sloperatio"]
        rltail = cfg["peakshape"]["lt_arearatio"]
        self.ltailslope_ratio = cfg["peakshape"]["lt_sloperatio"]
        rstep = cfg["peakshape"]["step_heightratio"]

        self.ratios = (rstail, rltail, rstep)

        setup.emin = cfg["fit"]["xmin"] * self.mcagain + self.mcazero
        setup.emax = cfg["fit"]["xmax"] * self.mcagain + self.mcazero


class sn3102(XRFDetector):
    aliases = ["XFlash5100"]

    def __init__(self, **kwargs):
        attenuators = kwargs.get("attenuators", {})
        ultralene = compoundfromname.compoundfromname("ultralene")
        moxtek = compoundfromname.compoundfromname("moxtek ap3.3")
        attenuators["FoilDetector"] = {
            "material": ultralene,
            "thickness": 4.064e-4,
        }  # cm
        attenuators["WindowDetector"] = {
            "material": moxtek,
            "thickness": 380e-4,
        }  # AP3.3 Moxtek
        attenuators["Detector"] = {
            "material": element.Element("Si"),
            "thickness": 500e-4,
        }
        kwargs["attenuators"] = attenuators

        activearea = kwargs.get("activearea", 0.8)
        kwargs["activearea"] = units.Quantity(activearea, "cm^2")
        kwargs["mcazero"] = kwargs.get("mcazero", 0.0)  # keV
        kwargs["mcagain"] = kwargs.get("mcagain", 5e-3)  # keV
        kwargs["mcanoise"] = kwargs.get("mcanoise", 0.1)  # keV
        kwargs["mcafano"] = kwargs.get("mcafano", 0.114)
        ehole = kwargs.get("ehole", constants.eholepair_si())
        kwargs["ehole"] = units.Quantity(ehole, "eV")
        super(sn3102, self).__init__(**kwargs)


class Leia(XRFDetector):
    aliases = ["SGX80"]

    def __init__(self, **kwargs):
        attenuators = kwargs.get("attenuators", {})
        ultralene = compoundfromname.compoundfromname("ultralene")
        attenuators["FoilDetector"] = {
            "material": ultralene,
            "thickness": 4.064e-4,
        }  # cm
        attenuators["WindowDetector"] = {
            "material": element.Element("Be"),
            "thickness": 25e-4,
        }
        attenuators["Detector"] = {
            "material": element.Element("Si"),
            "thickness": 450e-4,
        }
        kwargs["attenuators"] = attenuators

        activearea = kwargs.get("activearea", 0.8)
        kwargs["activearea"] = units.Quantity(activearea, "cm^2")
        kwargs["mcazero"] = kwargs.get("mcazero", 0.0)
        kwargs["mcagain"] = kwargs.get("mcagain", 5e-3)
        kwargs["mcanoise"] = kwargs.get("mcanoise", 50e-3)
        kwargs["mcafano"] = kwargs.get("mcafano", 0.19)
        ehole = kwargs.get("ehole", constants.eholepair_si())
        kwargs["ehole"] = units.Quantity(ehole, "eV")
        super(Leia, self).__init__(**kwargs)


class BB8(XRFDetector):
    aliases = ["SGX50"]

    def __init__(self, **kwargs):
        attenuators = kwargs.get("attenuators", {})
        ultralene = compoundfromname.compoundfromname("ultralene")
        attenuators["FoilDetector"] = {
            "material": ultralene,
            "thickness": 4.064e-4,
        }  # cm
        attenuators["WindowDetector"] = {
            "material": element.Element("Be"),
            "thickness": 12.5e-4,
        }
        attenuators["Detector"] = {
            "material": element.Element("Si"),
            "thickness": 450e-4,
        }
        kwargs["attenuators"] = attenuators

        activearea = kwargs.get("activearea", 0.5)
        kwargs["activearea"] = units.Quantity(activearea, "cm^2")
        kwargs["mcazero"] = kwargs.get("mcazero", 0.0)  # keV
        kwargs["mcagain"] = kwargs.get("mcagain", 5e-3)  # keV
        kwargs["mcanoise"] = kwargs.get("mcanoise", 0.1)  # keV
        kwargs["mcafano"] = kwargs.get("mcafano", 0.114)
        ehole = kwargs.get("ehole", constants.eholepair_si())
        kwargs["ehole"] = units.Quantity(ehole, "eV")
        super(BB8, self).__init__(**kwargs)


class DR40(XRFDetector):
    aliases = ["VITUSH80"]

    def __init__(self, **kwargs):
        attenuators = kwargs.get("attenuators", {})
        ultralene = compoundfromname.compoundfromname("ultralene")
        attenuators["FoilDetector"] = {
            "material": ultralene,
            "thickness": 4.064e-4,
        }  # cm
        attenuators["WindowDetector"] = {
            "material": element.Element("Be"),
            "thickness": 25e-4,
        }
        attenuators["Detector"] = {
            "material": element.Element("Si"),
            "thickness": 450e-4,
        }
        kwargs["attenuators"] = attenuators

        activearea = kwargs.get("activearea", 0.8)
        kwargs["activearea"] = units.Quantity(activearea, "cm^2")
        kwargs["mcazero"] = kwargs.get("mcazero", 0.0)  # keV
        kwargs["mcagain"] = kwargs.get("mcagain", 5e-3)  # keV
        kwargs["mcanoise"] = kwargs.get("mcanoise", 0.1)  # keV
        kwargs["mcafano"] = kwargs.get("mcafano", 0.114)
        ehole = kwargs.get("ehole", constants.eholepair_si())
        kwargs["ehole"] = units.Quantity(ehole, "eV")
        super(DR40, self).__init__(**kwargs)


class ID16b_Virtual1(XRFDetector):
    def __init__(self, **kwargs):
        attenuators = kwargs.get("attenuators", {})
        attenuators["WindowDetector"] = {
            "material": element.Element("Be"),
            "thickness": 50e-4,
        }
        attenuators["Deadlayer"] = {
            "material": element.Element("Si"),
            "thickness": 30e-7,
        }
        attenuators["Detector"] = {
            "material": element.Element("Si"),
            "thickness": 450e-4,
        }
        kwargs["attenuators"] = attenuators

        activearea = kwargs.get("activearea", 0.8)
        kwargs["activearea"] = units.Quantity(activearea, "cm^2")
        kwargs["mcazero"] = kwargs.get("mcazero", 0.0)  # keV
        kwargs["mcagain"] = kwargs.get("mcagain", 10e-3)  # keV
        kwargs["mcanoise"] = kwargs.get("mcagain", 0.1)  # keV
        kwargs["mcafano"] = kwargs.get("mcagain", 0.114)
        ehole = kwargs.get("ehole", constants.eholepair_si())
        kwargs["ehole"] = units.Quantity(ehole, "eV")
        super(ID16b_Virtual1, self).__init__(**kwargs)


factory = XRFDetector.factory
registry = XRFDetector.clsregistry

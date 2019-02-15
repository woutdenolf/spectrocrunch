# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from __future__ import division

from ..geometries import flatarea
from ..sources import xray as xraysources
from ..detectors import area
from ..materials import scintillators
from ..materials import lenses
from ..sources import emspectrum
from ..math import noisepropagation
from ..utils.instance import isarray

import numpy as np
import uncertainties
import itertools
import matplotlib.pyplot as plt


def transmission(N, N0, D=0, D0=0,
                 tframe_data=1, nframe_data=1,
                 tframe_flat=1, nframe_flat=1,
                 nframe_dark=1):
    num = N - D/nframe_dark*nframe_data
    num /= nframe_data*tframe_data
    denom = N0 - D0/nframe_dark*nframe_flat
    denom /= nframe_flat*tframe_flat
    return num/denom


def absorbance(_transmission):
    ret = np.empty_like(_transmission)
    ind = _transmission > 0
    ret[ind] = -uncertainties.unumpy.log(_transmission[ind])
    ret[~ind] = noisepropagation.randomvariable(-np.inf, np.inf)
    return ret


class id21_ffsetup(object):

    def __init__(self, sample=None, scint="LSO", lens="x10", **kwargs):
        self.setsample(sample, **kwargs)
        if scint == "LSO":
            self.oscint = scintillators.factory("LSO ID21", thickness=10)
        else:
            self.oscint = scintillators.factory("GGG ID21", thickness=13)
        if lens == "x10":
            self.olens = lenses.factory("mitutoyoid21_10x")
        else:
            self.olens = lenses.factory("mitutoyoid21_10x")

    def setsample(self, sample, det="pcoedge55"):
        self.sample = sample
        self.odetector = area.factory(det)

        if sample is not None:
            src = xraysources.factory("synchrotron")
            geometry = flatarea.factory(
                "perpendicular", detector=self.odetector, source=src)
            self.sample.geometry = geometry

    def propagate(self, N, E, tframe, nframe, samplein=False, withnoise=True, forward=True, poissonapprox=False):
        if isarray(E):
            E = np.asarray(E)
        bsample = samplein and self.sample is not None
        reshape = None
        if forward:
            if bsample:
                N = self.sample.propagate(N, E, forward=forward)
            N = self.oscint.propagate(N, E, forward=forward)
            if isarray(N):
                reshape = N.shape
                N = N.flatten()
            N = self.olens.propagate(N, self.oscint.visspectrum,
                                     nrefrac=self.oscint.get_nrefrac(),
                                     forward=forward)
            N = self.odetector.propagate(N, self.oscint.visspectrum,
                                         tframe=tframe, nframe=nframe,
                                         forward=forward, poissonapprox=poissonapprox)
            if reshape:
                N = N.reshape(reshape)
        else:
            if isarray(N):
                reshape = N.shape
                N = N.flatten()
            N = self.odetector.propagate(N, self.oscint.visspectrum,
                                         tframe=tframe, nframe=nframe,
                                         forward=forward)
            N = self.olens.propagate(N, self.oscint.visspectrum,
                                     nrefrac=self.oscint.get_nrefrac(),
                                     forward=forward)
            if reshape:
                N = N.reshape(reshape)

            N = self.oscint.propagate(N, E, forward=forward)
            if bsample:
                N = self.sample.propagate(N, E, forward=forward)
        return N

    def photontoDU(self, energy, tframe, nframe, samplein=False, withnoise=True):
        """Linear relation between incoming photons and DU

        Args:
            energy(num):

        Returns:
            num: DU/ph,DU
        """
        N0 = self.propagate(1, energy, tframe, nframe, samplein=samplein)
        D0 = self.odetector.propagate(
            0, emspectrum.Dark(), tframe=tframe, nframe=nframe)
        return N0-D0, D0

    def measurement(self, flux, energy,
                    tframe_data=None, nframe_data=None,
                    tframe_flat=None, nframe_flat=None,
                    nframe_dark=None, withnoise=True):
        """ID21 fullfield noise propagation

        Args:
            flux (num|array): incomming flux (ph/sec)
            energy (num|array): associated energy (keV)
            sample (spectrocrunch.simulation.materials.Material): sample
            tframe_data(num): time per frame (sec)
            nframe_data(num): number of data frames
            tframe_flat(num): time per frame (sec)
            nframe_flat(num): number of flat frames
            nframe_dark(num): number of dark frames
            scint(Optional(str)):
            lens(Optional(str)):

        Returns:
            4-tuple(uncertainties.unumpy.uarray): detector signal in DU (with and w.o. sample, with and w.o. exposure)
        """

        # With sample
        N = flux*tframe_data
        if withnoise:
            N = noisepropagation.poisson(N)
        N = self.propagate(N, energy, tframe_data, nframe_data, samplein=True)

        # Without sample
        N0 = flux*tframe_flat
        if withnoise:
            N0 = noisepropagation.poisson(N0)
        N0 = self.propagate(N0, energy, tframe_flat,
                            nframe_flat, samplein=False)

        # Without beam
        if withnoise:
            D = noisepropagation.randomvariable(0, 0)
            D0 = noisepropagation.randomvariable(0, 0)
        else:
            D = 0
            D0 = 0

        if nframe_dark != 0:
            D = self.odetector.propagate(
                D, self.oscint.visspectrum, tframe=tframe_data, nframe=nframe_dark)
            D0 = self.odetector.propagate(
                D0, self.oscint.visspectrum, tframe=tframe_flat, nframe=nframe_dark)

        return N, N0, D, D0

    def fluxfromflat(self, image, energy, tframe, nframe):
        return self.propagate(image, energy, tframe, nframe, samplein=False, withnoise=False, forward=False)/tframe

    def flatfromflux(self, flux, energy, tframe, nframe):
        return self.propagate(flux, energy, tframe, nframe, samplein=False, withnoise=False, forward=True)

    @staticmethod
    def getnframes(totaltime, tframe_data, tframe_flat, fracflat):
        ndata = max(int(round(totaltime/tframe_data*(1-fracflat))), 1)
        nflat = max(int(round(totaltime/tframe_data*fracflat/2.)), 1)
        nflat *= 2  # before and after
        return ndata, nflat

    def xanes(self, flux, energy, **kwarg):
        N, N0, D, D0 = self.measurement(flux, energy, **kwarg)
        T = transmission(N, N0, D=D, D0=D0, **kwarg)
        XAS = absorbance(T)
        signal = noisepropagation.E(XAS)
        noise = noisepropagation.S(XAS)
        return signal, noise

    def costfunc(self, flux, energy, **kwargs):
        signal, noise = self.xanes(flux, energy, **kwargs)
        jump = np.abs(signal[-1]-signal[0])
        return np.max(noise)/jump

    def __str__(self):
        return str(self.sample)

    def plotxanesnoise(self, flux, energy, **kwargs):
        signal, noise = self.xanes(flux, energy, **kwargs)
        signal = np.random.normal(signal, noise)
        p = plt.plot(energy, signal)
        plt.xlabel("Energy (keV)")
        plt.ylabel("Absorbance")
        return p

    def plotxanesNSR(self, flux, energy, **kwargs):
        signal, noise = self.xanes(flux, energy, **kwargs)
        p = plt.plot(energy, noise/signal*100)
        plt.xlabel("Energy (keV)")
        plt.ylabel("N/S (%)")
        return p

    def plotxanes(self, flux, energy, **kwargs):
        signal, _ = self.xanes(flux, energy, **kwargs)
        p = plt.plot(energy, signal)
        plt.xlabel("Energy (keV)")
        plt.ylabel("Absorbance")
        return p


def id21_transmissionnoise(flux, energy, time, iodetgain, idetgain, iodet="1"):
    """ID21 transmission noise propagation
    """
    pass


def id21_xrfnoise(flux, energy, time):
    """ID21 XRF noise propagation
    """
    pass

# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
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

from . import xrayspectrum
from ..utils import instance
from ..patch.pint import ureg
import numpy as np


def refractive_index_factor(energy, density):
    return ureg.Quantity(energy, 'keV').to("cm", "spectroscopy")**2 *\
          (ureg.classical_electron_radius*ureg.avogadro_number *
           ureg.Quantity(density, 'g/cm^3')/(2*np.pi))


def refractive_index_delta_calc(energy, e_wfrac, density, **kwargs):
    delta = sum(e_wfrac[e]/e.MM*e.scatfact_real(energy, **kwargs)
                for e in e_wfrac)
    delta = ureg.Quantity(delta, 'mol/g') *\
        refractive_index_factor(energy, density)
    return delta.to("dimensionless").magnitude


def refractive_index_beta_calc(energy, e_wfrac, density, **kwargs):
    # TODO: modify sign in scatfact_imag?
    beta = -sum(e_wfrac[e]/e.MM*e.scatfact_imag(energy, **kwargs)
                for e in e_wfrac)
    beta = ureg.Quantity(beta, 'mol/g') *\
        refractive_index_factor(energy, density)
    return beta.to("dimensionless").magnitude


class ElementBase(object):

    def refractive_index_delta(self, E, fine=False, decomposed=False, **kwargs):
        """n = 1-delta-i*beta
        """
        if hasattr(self, 'structure') and fine:
            environ = self
        else:
            environ = None
        return refractive_index_delta_calc(E, self.elemental_massfractions(),
                                           self.density, environ=environ, **kwargs)

    def refractive_index_beta(self, E, fine=False, decomposed=False, **kwargs):
        """n = 1-delta-i*beta
        """
        if hasattr(self, 'structure') and fine:
            environ = self
        else:
            environ = None
        return refractive_index_beta_calc(E, self.elemental_massfractions(),
                                          self.density, environ=environ, **kwargs)
        
    def refractive_index_real(self, E, **kwargs):
        """Real part of the refractive index
        """
        return 1-self.refractive_index_delta(E)

    def refractive_index_imag(self, E, **kwargs):
        """Imaginary part of the refractive index
        """
        return -self.refractive_index_beta(E)
        
    def xrayspectrum(self, E, source=None, weights=None, emin=0, emax=None, **kwargs):
        E = instance.asarray(E)
        if emax is None:
            emax = np.max(E)
        self.markabsorber(energybounds=(emin, emax))

        spectrum = xrayspectrum.Spectrum()

        if source is None:
            spectrum.update(self.fluorescence_cross_section_lines(
                E, decomposed=False, **kwargs))
            spectrum[xrayspectrum.RayleighLine(E)] = self.rayleigh_cross_section(
                E, decomposed=False, **kwargs)
            spectrum[xrayspectrum.ComptonLine(E)] = self.compton_cross_section(
                E, decomposed=False, **kwargs)
            spectrum.type = spectrum.TYPES.crosssection
        else:
            spectrum.density = self.density
            spectrum.update(self.diff_fluorescence_cross_section(
                E, decomposed=False, **kwargs))
            spectrum[xrayspectrum.RayleighLine(E)] = self.diff_rayleigh_cross_section(
                E, source=source, decomposed=False, **kwargs)
            spectrum[xrayspectrum.ComptonLine(E)] = self.diff_compton_cross_section(
                E, source=source, decomposed=False, **kwargs)
            spectrum.type = spectrum.TYPES.diffcrosssection

        spectrum.apply_weights(weights)
        spectrum.density = self.density
        spectrum.xlim = [emin, emax]
        spectrum.title = str(self)
        spectrum.geomkwargs = kwargs

        return spectrum

    def fisxgroups(self, emin=0, emax=np.inf):
        self.markabsorber(energybounds=[emin, emax])
        return {el: el.shells for el in self.elements}

    pymcamaterial_prefix = "Material_"
    pymcacomment_prefix = "From spectrocrunch: "

    @property
    def pymcaname(self):
        return self.pymcamaterial_prefix+self.name

    @property
    def pymcacomment(self):
        return self.pymcacomment_prefix+self.name

    @classmethod
    def namefrompymca(cls, string):
        for prefix in [cls.pymcamaterial_prefix, cls.pymcacomment_prefix]:
            if string.startswith(prefix):
                return string[len(prefix):]
        return string

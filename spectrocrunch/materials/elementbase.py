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
from ..common import instance

import numpy as np

class ElementBase(object):

    def xrayspectrum(self,E,source=None,weights=None,emin=0,emax=None,**kwargs):
        E = instance.asarray(E)
        if emax is None:
            emax = E[-1]
        self.markabsorber(energybounds=[emin,emax])
        
        spectrum = xrayspectrum.Spectrum()

        if source is None:
            spectrum.update(self.fluorescence_cross_section_lines(E,decomposed=False,**kwargs))
            spectrum[xrayspectrum.RayleighLine(E)] = self.rayleigh_cross_section(E,decomposed=False,**kwargs)
            spectrum[xrayspectrum.ComptonLine(E)] = self.compton_cross_section(E,decomposed=False,**kwargs)
            spectrum.type = spectrum.TYPES.crosssection
        else:
            spectrum.density = self.density
            spectrum.update(self.diff_fluorescence_cross_section(E,decomposed=False,**kwargs))
            spectrum[xrayspectrum.RayleighLine(E)] = self.diff_rayleigh_cross_section(E,source=source,decomposed=False,**kwargs)
            spectrum[xrayspectrum.ComptonLine(E)] = self.diff_compton_cross_section(E,source=source,decomposed=False,**kwargs)
            spectrum.type = spectrum.TYPES.diffcrosssection

        spectrum.apply_weights(weights)
        spectrum.density = self.density
        spectrum.xlim = [emin,emax]
        spectrum.title = str(self)
        
        return spectrum
        
    def fisxgroups(self,emin=0,emax=np.inf):
        self.markabsorber(energybounds=[emin,emax])
        return {el:el.shells for el in self.elements}
        

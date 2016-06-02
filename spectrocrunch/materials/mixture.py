# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
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

from .types import fractionType
from . import stoichiometry
import numpy as np

class mixture(object):

    def __init__(self,compounds,frac,fractype):
        """
        Args:
            compounds(list[obj]): list of compound objects
            frac(list[float]): compound fractions in the mixture
            fractype(fractionType): compound fraction type
        """
        
        # Mole fractions
        if fractype == fractionType.mole:
            MM = np.asarray([c.molarmass() for c in compounds])
            wfrac = stoichiometry.frac_mole_to_weight(np.asarray(frac),MM)
        elif fractype == fractionType.volume:
            rho = np.asarray([c.density for c in compounds])
            wfrac = stoichiometry.frac_volume_to_weight(np.asarray(frac),rho)
        else:
            wfrac = frac/sum(np.asarray(frac)) # sum is not necessarily zero in this case!

        self.compounds = dict(zip(compounds,wfrac))

    def __repr__(self):
        return '\n'.join("{} wt% {}".format(100*s[1],s[0]) for s in self.compounds.items())

    def __getitem__(self,compound):
        return self.compounds[compound]

    def density(self):
        MM = np.asarray([c.molarmass() for c in self.compounds])
        rho = np.asarray([c.density for c in self.compounds])
        wfrac = np.asarray(self.compounds.values())
        nfrac = stoichiometry.frac_weight_to_mole(wfrac,MM)
        return stoichiometry.density_from_molefrac(nfrac,rho,MM)

    def _crosssection(self,method,E,decimals=6,refresh=False,fine=False,decomposed=False):
        """Calculate compound cross-sections
        """
        if hasattr(self,'structure') and fine:
            environ = self
        else:
            environ = None

        # compound cross-sections
        if decomposed:
            ret = {}
            for c in self.compounds:
                if hasattr(c,'structure') and fine:
                    environ = c
                else:
                    environ = None

                wc = self.compounds[c]
                ret[c] = E*0
                for e in c.elements:
                    ret[c] += wc*c.elements[e]*getattr(e,method)(E,environ=environ,decimals=decimals,refresh=refresh)
        else:
            ret = E*0
            for c in self.compounds:
                if hasattr(c,'structure') and fine:
                    environ = c
                else:
                    environ = None

                wc = self.compounds[c]
                for e in c.elements:
                    ret += wc*c.elements[e]*getattr(e,method)(E,environ=environ,decimals=decimals,refresh=refresh)

        return ret

    def mass_att_coeff(self,E,decimals=6,refresh=False,fine=False,decomposed=False):
        """Mass attenuation coefficient (cm^2/g, E in keV). In other words: transmission XAS.
        """
        return self._crosssection("mass_att_coeff",E,decimals=decimals,refresh=refresh,fine=fine,decomposed=decomposed)

    def partial_mass_abs_coeff(self,E,decimals=6,refresh=False,fine=False,decomposed=False):
        """Mass absorption coefficient for the selected shells and lines (cm^2/g, E in keV). In other words: fluorescence XAS.
        """
        return self._crosssection("partial_mass_abs_coeff",E,decimals=decimals,refresh=refresh,fine=fine,decomposed=decomposed)

    def markabsorber(self,symb,shells=[],fluolines=[]):
        """
        Args:
            symb(str): element symbol
        """
        for c in self.compounds:
            c.markabsorber(symb,shells=shells,fluolines=fluolines)

    def unmarkabsorber(self):
        for c in self.compounds:
            c.markasabsorber()

    def get_energy(self,energyrange,defaultinc=1):
        """Get absolute energies (keV) from a relative energy range (eV)
        """
        for c in self.compounds:
            ret = c.get_energy(energyrange,defaultinc=defaultinc)
            if ret is not None:
                return ret
        return None



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
        
        # Compound mole fractions
        fractions = np.asarray(frac)
        if fractype == fractionType.mole:
            nfrac = fractions
        elif fractype == fractionType.volume:
            MM = np.asarray([c.molarmass() for c in compounds])
            rho = np.asarray([c.density for c in compounds])
            nfrac = stoichiometry.frac_volume_to_mole(fractions,rho,MM)
        else:
            MM = np.asarray([c.molarmass() for c in compounds])
            nfrac = stoichiometry.frac_weight_to_mole(fractions,MM)

        self.compounds = dict(zip(compounds,nfrac))

    def __repr__(self):
        return '\n'.join("{} {}".format(s[1],s[0]) for s in self.compounds.items())

    def __getitem__(self,compound):
        return self.compounds[compound]

    def molarmass(self):
        MM = np.asarray([c.molarmass() for c in self.compounds])
        nfrac = np.asarray(self.molefractions(total=True).values())
        return (MM*nfrac).sum()

    def weightfractions(self):
        MM = np.asarray([c.molarmass() for c in self.compounds])
        nfrac = np.asarray(self.molefractions().values())
        wfrac = stoichiometry.frac_mole_to_weight(nfrac,MM)
        return dict(zip(self.compounds.keys(),wfrac))

    def molefractions(self,total=True):
        if total:
            return self.compounds
        else:
            nfrac = np.asarray(self.compounds.values())
            nfrac /= nfrac.sum()
            return dict(zip(self.compounds.keys(),nfrac))

    def density(self):
        MM = np.asarray([c.molarmass() for c in self.compounds])
        rho = np.asarray([c.density for c in self.compounds])
        nfrac = np.asarray(self.molefractions().values())
        return stoichiometry.density_from_molefrac(nfrac,rho,MM)

    def elemental_molefractions(self,total=True):
        ret = {}

        c_nfrac = self.molefractions(total=total)
        for c in c_nfrac:
            e_nfrac = c.molefractions(total=total)
            for e in e_nfrac:
                if e in ret:
                    ret[e] += c_nfrac[c]*e_nfrac[e]
                else:
                    ret[e] = c_nfrac[c]*e_nfrac[e]
        return ret

    def elemental_weightfractions(self):
        ret = {}

        c_wfrac = self.weightfractions()
        for c in c_wfrac:
            e_wfrac = c.weightfractions()
            for e in e_wfrac:
                if e in ret:
                    ret[e] += c_wfrac[c]*e_wfrac[e]
                else:
                    ret[e] = c_wfrac[c]*e_wfrac[e]
        return ret

    def _crosssection(self,method,E,decimals=6,refresh=False,fine=False,decomposed=False):
        """Calculate compound cross-sections
        """
        if hasattr(self,'structure') and fine:
            environ = self
        else:
            environ = None

        # compound cross-sections
        c_wfrac = self.weightfractions()
        if decomposed:
            ret = {}
            for c in c_wfrac:
                if hasattr(c,'structure') and fine:
                    environ = c
                else:
                    environ = None

                e_wfrac = c.weightfractions()
                ret[c] = E*0
                for e in c.elements:
                    ret[c] += c_wfrac[c]*e_wfrac[e]*getattr(e,method)(E,environ=environ,decimals=decimals,refresh=refresh)
        else:
            ret = E*0
            for c in c_wfrac:
                if hasattr(c,'structure') and fine:
                    environ = c
                else:
                    environ = None

                e_wfrac = c.weightfractions()
                for e in c.elements:
                    ret += c_wfrac[c]*e_wfrac[e]*getattr(e,method)(E,environ=environ,decimals=decimals,refresh=refresh)

        return ret

    def mass_att_coeff(self,E,decimals=6,refresh=False,fine=False,decomposed=False):
        """Mass attenuation coefficient (cm^2/g, E in keV). In other words: transmission XAS.
        """
        return self._crosssection("mass_att_coeff",E,decimals=decimals,refresh=refresh,fine=fine,decomposed=decomposed)

    def mass_abs_coeff(self,E,decimals=6,refresh=False,fine=False,decomposed=False):
        """Mass absorption coefficient (cm^2/g, E in keV).
        """
        return self._crosssection("mass_abs_coeff",E,decimals=decimals,refresh=refresh,fine=fine,decomposed=decomposed)

    def partial_mass_abs_coeff(self,E,decimals=6,refresh=False,fine=False,decomposed=False):
        """Mass absorption coefficient for the selected shells and lines (cm^2/g, E in keV). In other words: fluorescence XAS.
        """
        return self._crosssection("partial_mass_abs_coeff",E,decimals=decimals,refresh=refresh,fine=fine,decomposed=decomposed)

    def scattering_cross_section(self,E,decimals=6,refresh=False,fine=False,decomposed=False):
        """Scattering cross section (cm^2/g, E in keV).
        """
        return self._crosssection("scattering_cross_section",E,decimals=decimals,refresh=refresh,fine=fine,decomposed=decomposed)

    def compton_cross_section(self,E,decimals=6,refresh=False,fine=False,decomposed=False):
        """Compton cross section (cm^2/g, E in keV).
        """
        return self._crosssection("compton_cross_section",E,decimals=decimals,refresh=refresh,fine=fine,decomposed=decomposed)

    def rayleigh_cross_section(self,E,decimals=6,refresh=False,fine=False,decomposed=False):
        """Rayleigh cross section (cm^2/g, E in keV).
        """
        return self._crosssection("rayleigh_cross_section",E,decimals=decimals,refresh=refresh,fine=fine,decomposed=decomposed)

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



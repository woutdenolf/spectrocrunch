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

from .compoundfromlist import compoundfromlist
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

        # Compounds (no duplicates)
        self._compose_compounds(compounds,nfrac)

    def _compose_compounds(self,compounds,nfrac):
        self.compounds = {}
        for c,n in zip(compounds,nfrac):
            if c in self.compounds:
                self.compounds[c] += float(n)
            else:
                self.compounds[c] = float(n)

    def addcompound(self,c,frac,fractype):
        """Add a compound to the mixture

        Args:
            c(compounds): compound
            frac(num): compound fraction
            fractype(fractionType): compound fraction type
        
        """
        if fractype == fractionType.mole:
            nfrac = np.asarray(self.molefractions().values())
            if nfrac.sum()==1:
                nfrac = stoichiometry.add_frac(nfrac,frac)
            else:
                nfrac = np.append(nfrac,frac)

            compounds = self.compounds.keys()+[c]
        elif fractype == fractionType.volume:
            vfrac = np.asarray(self.volumefractions().values())
            vfrac = stoichiometry.add_frac(vfrac,frac)

            compounds = self.compounds.keys()+[c]
            MM = np.asarray([c.molarmass() for c in compounds])
            rho = np.asarray([c.density for c in compounds])
            nfrac = stoichiometry.frac_volume_to_mole(vfrac,rho,MM)
        else:
            wfrac = np.asarray(self.weightfractions().values())
            wfrac = stoichiometry.add_frac(wfrac,frac)

            compounds = self.compounds.keys()+[c]
            MM = np.asarray([c.molarmass() for c in compounds])
            nfrac = stoichiometry.frac_weight_to_mole(wfrac,MM)

        # Compounds (no duplicates)
        self._compose_compounds(compounds,nfrac)

    def tocompound(self,name):
        tmp = self.elemental_molefractions()
        return compoundfromlist(tmp.keys(),tmp.values(),fractionType.mole,self.density,name=name)

    def __repr__(self):
        return '\n'.join("{} {}".format(s[1],s[0]) for s in self.compounds.items())

    def __getitem__(self,compound):
        return self.compounds[compound]

    def molarmass(self):
        MM = np.asarray([c.molarmass() for c in self.compounds])
        nfrac = np.asarray(self.molefractions(total=True).values())
        return (MM*nfrac).sum()

    def volumefractions(self):
        MM = np.asarray([c.molarmass() for c in self.compounds])
        rho = np.asarray([c.density for c in self.compounds])
        nfrac = np.asarray(self.molefractions().values())
        wfrac = stoichiometry.frac_mole_to_volume(nfrac,rho,MM)
        return dict(zip(self.compounds.keys(),wfrac))

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

    @property
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

    def _crosssection(self,method,E,fine=False,decomposed=False,**kwargs):
        """Calculate compound cross-sections
        """

        c_wfrac = self.weightfractions()
        if decomposed:
            ret = {}
            for c in c_wfrac:
                if hasattr(c,'structure') and fine:
                    environ = c
                else:
                    environ = None

                e_wfrac = c.weightfractions()

                ret[c] = {"frac":c_wfrac[c],"elements":{}}
                for e in c.elements:
                    ret[c]["elements"][e] += {"frac":e_wfrac[e],"cs":getattr(e,method)(E,environ=environ,**kwargs)}
        else:
            ret = E*0.
            for c in c_wfrac:
                if hasattr(c,'structure') and fine:
                    environ = c
                else:
                    environ = None

                e_wfrac = c.weightfractions()
                for e in c.elements:
                    ret += c_wfrac[c]*e_wfrac[e]*getattr(e,method)(E,environ=environ,**kwargs)

        return ret

    def mass_att_coeff(self,E,fine=False,decomposed=False,**kwargs):
        """Mass attenuation coefficient (cm^2/g, E in keV). Use for transmission XAS.
        """
        return self._crosssection("mass_att_coeff",E,fine=fine,decomposed=decomposed,**kwargs)

    def mass_abs_coeff(self,E,fine=False,decomposed=False,**kwargs):
        """Mass absorption coefficient (cm^2/g, E in keV).
        """
        return self._crosssection("mass_abs_coeff",E,fine=fine,decomposed=decomposed,**kwargs)

    def partial_mass_abs_coeff(self,E,fine=False,decomposed=False,**kwargs):
        """Mass absorption coefficient for the selected shells and lines (cm^2/g, E in keV).
        """
        return self._crosssection("partial_mass_abs_coeff",E,fine=fine,decomposed=decomposed,**kwargs)

    def scattering_cross_section(self,E,fine=False,decomposed=False,**kwargs):
        """Scattering cross section (cm^2/g, E in keV).
        """
        return self._crosssection("scattering_cross_section",E,fine=fine,decomposed=decomposed,**kwargs)

    def compton_cross_section(self,E,fine=False,decomposed=False,**kwargs):
        """Compton cross section (cm^2/g, E in keV).
        """
        return self._crosssection("compton_cross_section",E,fine=fine,decomposed=decomposed,**kwargs)

    def rayleigh_cross_section(self,E,fine=False,decomposed=False,**kwargs):
        """Rayleigh cross section (cm^2/g, E in keV).
        """
        return self._crosssection("rayleigh_cross_section",E,fine=fine,decomposed=decomposed,**kwargs)

    def xrf_cross_section(self,E,fine=False,decomposed=False,**kwargs):
        """XRF cross section (cm^2/g, E in keV). Use for fluorescence XAS.
        """
        return self._crosssection("xrf_cross_section",E,fine=fine,decomposed=decomposed,**kwargs)

    def xrf_cross_section_decomposed(self,E,fine=False,**kwargs):
        """XRF cross section (cm^2/g, E in keV).
        """
        return self._crosssection("xrf_cross_section_decomposed",E,fine=fine,decomposed=True,**kwargs)

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



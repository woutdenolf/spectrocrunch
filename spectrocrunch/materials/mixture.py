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

import numpy as np
import re

from . import compound
from . import compoundfromformula
from . import compoundfromlist
from . import element
from .types import fractionType
from . import stoichiometry
from . import xrayspectrum
from ..common import instance

class Mixture(object):

    def __init__(self,compounds,frac,fractype,name=None):
        """
        Args:
            compounds(list[obj]): list of compound objects
            frac(list[float]): compound fractions in the mixture
            fractype(fractionType): compound fraction type
            name(Optional[str]): compound name
        """
        
        # Mixture name
        if name is None:
            self.name = str(id(self))
        else:
            self.name = name
        
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
        
    def change_fractions(self,dfrac,fractype):
        """Change the compound fractions

        Args:
            frac(dict): compound fractions
            fractype(fractionType): element fraction type
        """
        compounds = dfrac.keys()
        fractions = np.asarray(dfrac.values())
        
        if fractype == fractionType.mole:
            nfrac = fractions
        elif fractype == fractionType.volume:
            MM = np.asarray([c.molarmass() for c in compounds])
            rho = np.asarray([c.density for c in compounds])
            nfrac = stoichiometry.frac_volume_to_mole(fractions,rho,MM)
        else:
            MM = np.asarray([c.molarmass() for c in compounds])
            nfrac = stoichiometry.frac_weight_to_mole(fractions,MM)
            
        self._compose_compounds(compounds,nfrac)
        
    def tocompound(self,name):
        tmp = self.elemental_molefractions()
        return compoundfromlist.CompoundFromList(tmp.keys(),tmp.values(),fractionType.mole,self.density,name=name)

    def __str__(self):
        ws = self.weightfractions()
        return ' + '.join("{:.02f} wt% {}".format(s[1]*100,s[0]) for s in ws.items())

    def __getitem__(self,compound):
        return self.compounds[compound]

    @property
    def elements(self):
        return list(set(e for c in self.compounds for e in c.elements))

    def molarmass(self,total=True):
        # equivalent to using molefractions when total=True
        # not equivalent to using molefractions when total=False
        nfrac = self.elemental_molefractions(total=total)
        MM = np.asarray([c.molarmass() for c in nfrac])
        nfrac = np.asarray(nfrac.values())
        return (MM*nfrac).sum()

    def molarmasseff(self):
        return self.molarmass(total=False)

    @property
    def Zeff(self):
        # equivalent to using elemental_molefractions
        nfrac = self.molefractions(total=False)
        Z = np.asarray([c.Zeff for c in nfrac])
        nfrac = np.asarray(nfrac.values())
        return (Z*nfrac).sum()
    
    def volumefractions(self):
        nfrac = self.molefractions()
        MM = np.asarray([c.molarmass() for c in nfrac])
        rho = np.asarray([c.density for c in nfrac])
        wfrac = stoichiometry.frac_mole_to_volume(np.asarray(nfrac.values()),rho,MM)
        return dict(zip(nfrac.keys(),wfrac))

    def weightfractions(self):
        nfrac = self.molefractions()
        MM = np.asarray([c.molarmass() for c in nfrac])
        wfrac = stoichiometry.frac_mole_to_weight(np.asarray(nfrac.values()),MM)
        return dict(zip(nfrac.keys(),wfrac))

    def molefractions(self,total=True):
        if total:
            return dict(self.compounds)
        else:
            nfrac = np.asarray(self.compounds.values())
            nfrac /= nfrac.sum()
            return dict(zip(self.compounds.keys(),nfrac))

    @property
    def nelements(self):
        return sum(c.nelements for c in self.compounds)
    
    @property
    def ncompounds(self):
        return len(self.compounds)
        
    def arealdensity(self):
        """Areal density in ng/mm^2
        """
        arealdensity = {}
        for c in self.compounds:
            tmp = c.arealdensity()
            for e,ad in tmp.items():
                if e in arealdensity:
                    arealdensity[e] += ad
                else:
                    arealdensity[e] = ad
        
        return arealdensity
        
    @property
    def density(self):
        nfrac = self.molefractions()
        MM = np.asarray([c.molarmass() for c in nfrac])
        rho = np.asarray([c.density for c in nfrac])
        nfrac = np.asarray(nfrac.values())
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

    @staticmethod
    def _cs_scattering(method):
        return method=="scattering_cross_section" or method=="compton_cross_section" or method=="rayleigh_cross_section"
    
    @staticmethod
    def _cs_dict(method):
        return method=="fluorescence_cross_section_lines"
        
    def _crosssection(self,method,E,fine=False,decomposed=False,**kwargs):
        """Calculate compound cross-sections
        """
        
        bscat = self._cs_scattering(method)
        if bscat and not self.hasscatterers:
            if decomposed:
                return {}
            else:
                return E*0.

        c_wfrac = self.weightfractions()
        if decomposed:
            ret = {}
            for c in c_wfrac:
                if bscat and not c.isscatterer:
                    continue
                
                ret[c] = {"w":c_wfrac[c],"cs":{}}
                
                if hasattr(c,'structure') and fine:
                    environ = c
                else:
                    environ = None

                e_wfrac = c.weightfractions()
                for e in e_wfrac:
                    ret[c]["cs"][e] = {"w":e_wfrac[e],"cs":getattr(e,method)(E,environ=environ,**kwargs)}
        else:
            if self._cs_dict(method):
                ret = {}
                
                for c in c_wfrac:
                    if self._cs_scattering(method) and not c.isscatterer:
                        continue
                        
                    if hasattr(c,'structure') and fine:
                        environ = c
                    else:
                        environ = None

                    e_wfrac = c.weightfractions()
                    for e in c.elements:
                        cs = getattr(e,method)(E,environ=environ,**kwargs)
                        if not cs:
                            continue
                            
                        w = c_wfrac[c]*e_wfrac[e]
                        for k,v in cs.items():
                            if k in ret:
                                ret[k] += w*v
                            else:
                                ret[k] = w*v
            else:
                ret = E*0.
                
                for c in c_wfrac:
                    if self._cs_scattering(method) and not c.isscatterer:
                        continue
                        
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

    def fluorescence_cross_section(self,E,fine=False,decomposed=False,**kwargs):
        """XRF cross section (cm^2/g, E in keV). Use for fluorescence XAS.
        """
        return self._crosssection("fluorescence_cross_section",E,fine=fine,decomposed=decomposed,**kwargs)

    def fluorescence_cross_section_lines(self,E,fine=False,decomposed=False,**kwargs):
        """XRF cross section (cm^2/g, E in keV).  Use for XRF.
        """
        return self._crosssection("fluorescence_cross_section_lines",E,fine=fine,decomposed=decomposed,**kwargs)
    
    def xrayspectrum(self,E,weights=None,emin=0,emax=None):
        E = instance.asarray(E)
        if emax is None:
            emax = E[-1]
        self.markabsorber(energybounds=[emin,emax])
        
        spectrum = xrayspectrum.Spectrum()
        spectrum.density = self.density
        spectrum.update(self.fluorescence_cross_section_lines(E,decomposed=False))
        spectrum[xrayspectrum.RayleighLine(E)] = self.rayleigh_cross_section(E,decomposed=False)
        spectrum[xrayspectrum.ComptonLine(E)] = self.compton_cross_section(E,decomposed=False)
        spectrum.apply_weights(weights)
        spectrum.xlim = [emin,emax]
        spectrum.title = str(self)
        spectrum.type = spectrum.TYPES.crosssection

        return spectrum
    
    def refractive_index_re(self,E,**kwargs):
        return 1-self.refractive_index_delta(E)
        
    def refractive_index_im(self,E,**kwargs):
        return -self.refractive_index_beta(E)
        
    def refractive_index_delta(self,E,fine=False,decomposed=False,**kwargs):
        return compound.Compound.refractive_index_delta_calc(E,self.elemental_weightfractions(),self.density,environ=None,**kwargs)
    
    def refractive_index_beta(self,E,fine=False,decomposed=False,**kwargs):
        return compound.Compound.refractive_index_beta_calc(E,self.elemental_weightfractions(),self.density,environ=None,**kwargs)

    def markinfo(self):
        yield "{}".format(self.name)
        for c in self.compounds:
            for s in c.markinfo():
                yield " {}".format(s)
        
    def markabsorber(self,symb=None,shells=None,fluolines=None,energybounds=None):
        """
        Args:
            symb(str): element symbol
        """
        for c in self.compounds:
            c.markabsorber(symb,shells=shells,fluolines=fluolines,energybounds=energybounds)

    def unmarkabsorber(self):
        for c in self.compounds:
            c.unmarkabsorber()

    def hasabsorbers(self):
        return any([c.hasabsorbers() for c in self.compounds])
        
    def markscatterer(self,name=None):
        """
        Args:
            name(str): compound name
        """
        for c in self.compounds:
            c.markscatterer(name)

    def unmarkscatterer(self):
        for c in self.compounds:
            c.unmarkscatterer()

    def hasscatterers(self):
        return any([c.isscatterer() for c in self.compounds])
            
    def get_energy(self,energyrange,defaultinc=1):
        """Get absolute energies (keV) from a relative energy range (eV)
        """
        for c in self.compounds:
            ret = c.get_energy(energyrange,defaultinc=defaultinc)
            if ret is not None:
                return ret
        return None

    @property
    def pymcaname(self):
        return self.name

    def topymca(self,defaultthickness=1e-4):
        return self.tocompound(self.pymcaname).topymca(defaultthickness=defaultthickness)
    
    def tofisx(self):
        return self.tocompound(self.pymcaname).tofisx()
        
    @classmethod
    def frompymca(cls,dic):
        # dic will always be a compound or mixture with a name, not a chemical formula
        
        # Assume all individual compounds have the density of the mixture
        purelement = "^(?P<element>[A-Z][a-z]?)1$"

        if instance.isstring(dic["CompoundList"]):
            lst = [dic["CompoundList"]]
        else:
            lst = dic["CompoundList"]
            
        n = len(lst)
        elements = [""]*n
        for i,c in enumerate(lst):
            m = re.match(purelement,c)
            if m:
                m = m.groupdict()
                elements[i] = m["element"]

        if all(elements):
            return compoundfromlist.CompoundFromList(elements,dic["CompoundFraction"],fractionType.weight,dic["Density"],name=dic["Comment"])
        else:
            compounds = [compoundfromformula.CompoundFromFormula(c,dic["Density"]) for c in lst]
            wfrac = dic["CompoundFraction"]
            return cls(compounds,wfrac,fractionType.weight,name=dic["Comment"])

    def fisxgroups(self,emin=0,emax=np.inf):
        self.markabsorber(energybounds=[emin,emax])
        return {el:el.shells for el in self.elements}


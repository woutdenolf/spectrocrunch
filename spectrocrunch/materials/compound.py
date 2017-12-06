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

from . import element
from . import interaction
from .types import fractionType
from . import stoichiometry

from ..common.hashable import Hashable
from ..common import listtools
from .. import ureg

import numpy as np

class Compound(Hashable):
    """Interface to a compound
    """

    def __init__(self,elements,frac,fractype,density,nrefrac=1,name=None):
        """
        Args:
            elements(list): list of elements (["Fe","O"] or [element("Fe"),element("O")])
            frac(list[float]): element fractions
            fractype(fractionType): element fraction type
            density(num): compound density in g/cm^3
            nrefrac(num): refractive index
            name(Optional[str]): compound name
        """

        # Compound name
        if name is None:
            self.name = str(id(self))
        else:
            self.name = name
        self.nrefrac = float(nrefrac)

        # Element mole fractions
        if fractype == fractionType.mole:
            nfrac = frac # keep unnormalized!
        elif fractype == fractionType.volume:
            # would be possible if you give the element densities in the compound
            # (which is not the same as the pure element density) but that's an unlikely given
            raise ValueError("Cannot create a compound from elemental volume fractions")
        else:
            elements = [element.Element(e) for e in elements]
            MM = np.asarray([e.MM for e in elements])
            nfrac = stoichiometry.frac_weight_to_mole(np.asarray(frac),MM) # normalized

        # Elements (no duplicates)
        self._compose_elements(elements,nfrac)

        # Compound density
        self.density = float(density)
        if self.density==0:
            if len(self.elements)==1:
                self.density = self.elements.keys()[0].density
            else:
                #rho = [e.density for e in self.elements]
                #self.density = np.mean(rho) # not based on anything, just a value
                if len(self.elements)==0:
                    self.density = 0.
                else:
                    self.density = 1. # approx. density of water

        self.isscatterer = True

    def _compose_elements(self,elements,nfrac):
        self.elements = {}
        for e,n in zip(elements,nfrac):
            if not isinstance(e,element.Element):
                e = element.Element(e)
            if e in self.elements:
                self.elements[e] += float(n)
            else:
                self.elements[e] = float(n)

    def addelement(self,el,frac,fractype,density=None):
        """Add an element to the compound

        Args:
            elements(str|element): "Fe" or element("Fe")
            frac(num): element fraction
            fractype(fractionType): element fraction type
            density(Optional(num)): new compound density in g/cm^3
        """

        if fractype == fractionType.mole:
            nfrac = np.asarray(self.molefractions().values())
            if nfrac.sum()==1:
                nfrac = stoichiometry.add_frac(nfrac,frac)
            else:
                nfrac = np.append(nfrac,frac)

            elements = self.elements.keys()+[element.Element(el)]
        elif fractype == fractionType.volume:
            raise ValueError("Cannot create a compound from elemental volume fractions")
        else:
            wfrac = np.asarray(self.weightfractions().values())
            wfrac = stoichiometry.add_frac(wfrac,frac)

            elements = self.elements.keys()+[element.Element(el)]
            MM = np.asarray([e.MM for e in elements])
            nfrac = stoichiometry.frac_weight_to_mole(wfrac,MM) # normalized

        # Elements (no duplicates)
        self._compose_elements(elements,nfrac)

    def change_fractions(self,frac,fractype):
        """Change the element fractions

        Args:
            frac(array): element fractions
            fractype(fractionType): element fraction type
        """
    
        elements = self.elements.keys()
    
        # Element mole fractions
        if fractype == fractionType.mole:
            nfrac = frac # keep unnormalized!
        elif fractype == fractionType.volume:
            # would be possible if you give the element densities in the compound
            # (which is not the same as the pure element density) but that's an unlikely given
            raise ValueError("Cannot create a compound from elemental volume fractions")
        else:
            MM = np.asarray([e.MM for e in elements])
            nfrac = stoichiometry.frac_weight_to_mole(np.asarray(frac),MM) # normalized
            
        # Elements
        self._compose_elements(elements,nfrac)

    def _cmpkey(self):
        return self.name

    def _stringrepr(self):
        return self.name
        return '{}: {} ({} g/cm^3)'.format(\
                    self.name,\
                    ', '.join("{} {}".format(s[1],s[0]) for s in sorted(self.elements.items())),\
                    self.density )

    def __getitem__(self,el):
        return self.elements[el]

    def molarmass(self):
        MM = np.asarray([e.MM for e in self.elements])
        nfrac = np.asarray(self.molefractions(total=True).values())
        return (MM*nfrac).sum()

    def weightfractions(self):
        MM = np.asarray([e.MM for e in self.elements])
        nfrac = np.asarray(self.molefractions().values())
        wfrac = stoichiometry.frac_mole_to_weight(nfrac,MM)
        return dict(zip(self.elements.keys(),wfrac))

    def molefractions(self,total=True):
        if total:
            return dict(self.elements)
        else:
            nfrac = np.asarray(self.elements.values())
            nfrac /= nfrac.sum()
            return dict(zip(self.elements.keys(),nfrac))

    @property
    def nelements(self):
        return len(self.elements)

    def arealdensity(self):
        """Areal density in ng/mm^2
        """
        # arealdensity (ng/mm^2) = w * density(g/cm^3) * 2e-5 (cm) * 1e-2 (cm^2/mm^2) * 1e9 (ng/g)
        w = self.weightfractions()
        arealdensity = np.asarray(w.values())*(self.density*200)
        return dict(zip(w.keys(),arealdensity))
        
    def markabsorber(self,symb,shells=[],fluolines=[]):
        """
        Args:
            symb(str): element symbol
        """
        for e in self.elements:
            e.markabsorber(symb,shells=shells,fluolines=fluolines)

    def unmarkabsorber(self):
        for e in self.elements:
            e.unmarkabsorber()

    def hasabsorbers(self):
        return any([e.isabsorber() for e in self.elements])

    def markscatterer(self,name):
        """
        Args:
            name(str): compound name
        """
        self.isscatterer = self==name

    def unmarkscatterer(self):
        self.isscatterer = False
        
    def getabsorberinfo(self):
        ret = {}
        for e in self.elements:
            if e.isabsorber() and e not in ret:
                ret[e] = e.getfluoinfo()
        return ret

    @staticmethod
    def _cs_scattering(method):
        return method=="scattering_cross_section" or method=="compton_cross_section" or method=="rayleigh_cross_section"

    def _crosssection(self,method,E,fine=False,decomposed=False,**kwargs):
        """Calculate compound cross-sections
        """
        
        if self._cs_scattering(method) and not self.isscatterer:
            if decomposed:
                return {}
            else:
                return E*0.

        if hasattr(self,'structure') and fine:
            environ = self
        else:
            environ = None

        e_wfrac = self.weightfractions()
        if decomposed:
            ret = {}
            for e in e_wfrac:
                ret[e] = {"w":e_wfrac[e],"cs":getattr(e,method)(E,environ=environ,**kwargs)}
        else:
            ret = E*0.
            for e in e_wfrac:
                ret += e_wfrac[e]*getattr(e,method)(E,environ=environ,**kwargs)

        return ret

    def mass_att_coeff(self,E,fine=False,decomposed=False,**kwargs):
        """Mass attenuation coefficient (cm^2/g, E in keV). Use for transmission XAS.
        """
        return self._crosssection("mass_att_coeff",E,fine=fine,decomposed=decomposed,**kwargs)

    def mass_abs_coeff(self,E,fine=False,decomposed=False,**kwargs):
        """Mass absorption coefficient (cm^2/g, E in keV)
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

    def refractive_index_re(self,E,fine=False,decomposed=False,**kwargs):
        """
        """
        if hasattr(self,'structure') and fine:
            environ = self
        else:
            environ = None
        
        e_wfrac = self.weightfractions()
        
        ret = E*0.
        
        for e in e_wfrac:
            ret += e_wfrac[e]*e.scatfact_re(E,environ=environ,**kwargs)/e.MM
        ret = 1 - ureg.Quantity(ret,'mol/g') *\
                  ureg.Quantity(E,'keV').to("cm","spectroscopy")**2 *\
                  (ureg.re*ureg.avogadro_number*ureg.Quantity(self.density,'g/cm^3')/(2*np.pi))
        
        return ret.to("dimensionless").magnitude

    def refractive_index_im(self,E,fine=False,decomposed=False,**kwargs):
        """
        """
        if hasattr(self,'structure') and fine:
            environ = self
        else:
            environ = None
        
        e_wfrac = self.weightfractions()
        
        ret = E*0.
        
        for e in e_wfrac:
            ret += e_wfrac[e]*e.scatfact_im(E,environ=environ,**kwargs)/e.MM
        ret = ureg.Quantity(ret,'mol/g') *\
              ureg.Quantity(E,'keV').to("cm","spectroscopy")**2 *\
              (-ureg.re*ureg.avogadro_number*ureg.Quantity(self.density,'g/cm^3')/(2*np.pi))
        
        return ret.to("dimensionless").magnitude
        
    def get_energy(self,energyrange,defaultinc=1):
        """Get absolute energies (keV) from a relative energy range (eV)

        Args:
            energyrange(np.array): energy range relative to the absorber's edges (if any)
        """
        for e in self.elements:
            ret = e.get_energy(energyrange,defaultinc=defaultinc)
            if ret is not None:
                return ret
        return None

    def pymcaformat(self):
        r = self.weightfractions()
        value = {'Comment': self.name,
                'CompoundFraction': r.values(),
                'Thickness': 0.,
                'Density': self.density,
                'CompoundList': ['{}1'.format(e) for e in r]}
        return self.name,value

    def fluointeractions(self):
        """Fluorescence interactions
        """
        for e in self.elements:
            for s in e.shells:
                for l in s.fluolines:
                    fl = interaction.InteractionFluo(e,s,l)
                    if fl.energy!=0:
                        yield fl
                        

        

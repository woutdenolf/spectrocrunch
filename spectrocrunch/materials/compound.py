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
import fisx

from . import element
from . import interaction
from .types import fractionType
from . import stoichiometry
from ..common.hashable import Hashable
from ..common import listtools
from ..common import instance
from .. import ureg
from . import xrayspectrum


class Compound(Hashable):
    """Interface to a compound
    """

    def __init__(self,elements,frac,fractype,density=None,nrefrac=1,name=None):
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
        nfrac = self.molefractions(total=True)
        MM = np.asarray([e.MM for e in nfrac])
        nfrac = np.asarray(nfrac.values())
        return (MM*nfrac).sum()

    def weightfractions(self):
        nfrac = self.molefractions()
        MM = np.asarray([e.MM for e in nfrac])
        wfrac = stoichiometry.frac_mole_to_weight(np.asarray(nfrac.values()),MM)
        return dict(zip(nfrac.keys(),wfrac))

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

    @property
    def ncompounds(self):
        return 1
        
    def arealdensity(self):
        """Areal density in ng/mm^2
        """
        # arealdensity (ng/mm^2) = w * density(g/cm^3) * 2e-5 (cm) * 1e-2 (cm^2/mm^2) * 1e9 (ng/g)
        w = self.weightfractions()
        arealdensity = np.asarray(w.values())*(self.density*200)
        return dict(zip(w.keys(),arealdensity))
        
    def markabsorber(self,symb=None,shells=None,fluolines=None,energybounds=None):
        """
        Args:
            symb(str): element symbol
        """
        for e in self.elements:
            if energybounds is None:
                energybounds2 = None
            else:
                energybounds2 = energybounds
            e.markabsorber(symb=symb,shells=shells,fluolines=fluolines,energybounds=energybounds2)

    def unmarkabsorber(self):
        for e in self.elements:
            e.unmarkabsorber()

    def hasabsorbers(self):
        return any([e.isabsorber() for e in self.elements])

    def markscatterer(self,name=None):
        """
        Args:
            name(str): compound name
        """
        if name is None:
            self.isscatterer = True
        else:
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

    @staticmethod
    def _cs_dict(method):
        return method=="fluorescence_cross_section_lines"
        
    def _crosssection(self,method,E,fine=False,decomposed=False,**kwargs):
        """Calculate compound cross-sections
        """
        
        if self._cs_scattering(method) and not self.isscatterer:
            if decomposed:
                return {}
            else:
                return np.zeros_like(E,dtype=float)

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
            if self._cs_dict(method):
                ret = {}
                for e in e_wfrac:
                    cs = getattr(e,method)(E,environ=environ,**kwargs)
                    if not cs:
                        continue
                    w = e_wfrac[e]
                    for k,v in cs.items():
                        if k in ret:
                            ret[k] += w*v
                        else:
                            ret[k] = w*v
            else:
                ret = np.zeros_like(E,dtype=float)
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

    def fluorescence_cross_section_lines(self,E,fine=False,decomposed=False,**kwargs):
        """XRF cross section (cm^2/g, E in keV). Use for XRF.
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

    @property
    def pymcaname(self):
        return self.name

    def topymca(self,defaultthickness=1e-4):
        r = self.weightfractions()
        value = {'Comment': self.pymcaname,
                'CompoundFraction': r.values(),
                'Thickness': defaultthickness,
                'Density': self.density,
                'CompoundList': ['{}1'.format(e) for e in r]}
        return self.pymcaname,value

    def tofisx(self):
        r = self.weightfractions()
        o = fisx.Material(self.pymcaname, self.density, 1e-10)
        o.setCompositionFromLists(['{}1'.format(e) for e in r],r.values())
        return o
        
    def fisxgroups(self,emin=0,emax=np.inf):
        self.markabsorber(energybounds=[emin,emax])
        return {el:el.shells for el in self.elements}
        

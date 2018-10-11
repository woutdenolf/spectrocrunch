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

from . import multielementbase
from . import element
from . import interaction
from . import types
from . import stoichiometry
from ..utils.hashable import Hashable
from ..utils import listtools
from ..utils import instance

import numpy as np
import fisx

class Compound(Hashable,multielementbase.MultiElementBase):
    """Interface to a compound
    """

    def __init__(self,elements,frac,fractype,density=None,nrefrac=1,name=None):
        """
        Args:
            elements(list): list of elements (["Fe","O"] or [element("Fe"),element("O")])
            frac(list[float]): element fractions
            fractype(types.fraction): element fraction type
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
        if fractype == types.fraction.mole:
            nfrac = frac # keep unnormalized!
        elif fractype == types.fraction.volume:
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
        if density==0 or density is None:
            if len(self._elements)==1:
                self.density = self._elements.keys()[0].density
            else:
                #rho = [e.density for e in self._elements]
                #self.density = np.mean(rho) # not based on anything, just a value
                if len(self._elements)==0:
                    self.density = 0.
                else:
                    self.density = 1. # approx. density of water
        else:
            self.density = float(density)
            
        self.isscatterer = True

    def _compose_elements(self,elements,nfrac):
        self._elements = {}
        for e,n in zip(elements,nfrac):
            if not isinstance(e,element.Element):
                e = element.Element(e)
            if e in self._elements:
                self._elements[e] += float(n)
            else:
                self._elements[e] = float(n)

    def addelements(self,elements,fractions,fractype,density=None):
        """Add an element to the compound

        Args:
            elements(str|element): "Fe" or element("Fe")
            frac(num): element fraction
            fractype(types.fraction): element fraction type
            density(Optional(num)): new compound density in g/cm^3
        """

        if fractype == types.fraction.volume:
            raise ValueError("Cannot create a compound from elemental volume fractions")
            
        if not instance.isarray(elements):
            elements = [elements]
        if not instance.isarray(fractions):
            fractions = [fractions]
            
        if fractype == types.fraction.mole:
            nfrac = np.asarray(self.molefractions().values())
            if nfrac.sum()==1:
                nfrac = stoichiometry.add_frac(nfrac,fractions)
            else:
                nfrac = np.append(nfrac,fractions)

            elements = self._elements.keys()+[element.Element(el) for el in elements]
        else:
            wfrac = np.asarray(self.massfractions().values())
            wfrac = stoichiometry.add_frac(wfrac,fractions)

            elements = self._elements.keys()+[element.Element(el) for el in elements]
            MM = np.asarray([e.MM for e in elements])
            nfrac = stoichiometry.frac_weight_to_mole(wfrac,MM) # normalized

        # Elements (no duplicates)
        self._compose_elements(elements,nfrac)

    def change_fractions(self,dfrac,fractype):
        """Change the element fractions

        Args:
            dfrac(dict): element fractions
            fractype(types.fraction): fraction type
        """
    
        elements = dfrac.keys()
        fractions = np.asarray(dfrac.values())
        
        # Element mole fractions
        if fractype == types.fraction.mole:
            nfrac = fractions # keep unnormalized!
        elif fractype == types.fraction.volume:
            # would be possible if you give the element densities in the compound
            # (which is not the same as the pure element density) but that's an unlikely given
            raise ValueError("Cannot create a compound from elemental volume fractions")
        else:
            MM = np.asarray([e.MM for e in elements])
            nfrac = stoichiometry.frac_weight_to_mole(fractions,MM) # normalized
            
        # Elements
        self._compose_elements(elements,nfrac)

    def _cmpkey(self):
        return self.name

    def _stringrepr(self):
        return self.name
        return '{}: {} ({} g/cm^3)'.format(\
                    self.name,\
                    ', '.join("{} {}".format(s[1],s[0]) for s in sorted(self._elements.items())),\
                    self.density )

    def __getitem__(self,el):
        return self._elements[el]

    @property
    def elements(self):
        return self._elements.keys()

    @property
    def parts(self):
        return self._elements
        
    def molarmass(self,total=True):
        nfrac = self.molefractions(total=total)
        MM = np.asarray([e.MM for e in nfrac])
        nfrac = np.asarray(nfrac.values())
        return (MM*nfrac).sum()

    def molarmasseff(self):
        return self.molarmass(total=False)

    @property
    def Zeff(self):
        nfrac = self.molefractions(total=False)
        Z = np.asarray([e.Z for e in nfrac])
        nfrac = np.asarray(nfrac.values())
        return (Z*nfrac).sum()

    def molefractions(self,total=True):
        if total:
            return dict(self._elements)
        else:
            nfrac = np.asarray(self._elements.values())
            nfrac /= nfrac.sum()
            return dict(zip(self._elements.keys(),nfrac))
            
    def massfractions(self):
        nfrac = self.molefractions()
        MM = np.asarray([e.MM for e in nfrac])
        wfrac = stoichiometry.frac_mole_to_weight(np.asarray(nfrac.values()),MM)
        return dict(zip(nfrac.keys(),wfrac))

    def elemental_molefractions(self,**kwargs):
        return self.molefractions(**kwargs)
        
    def elemental_massfractions(self):
        return self.massfractions()

    def fractions(self,fractype):
        if fractype == types.fraction.mole:
            return self.molefractions()
        elif fractype == types.fraction.volume:
            raise ValueError("Cannot calculate elemental volume fractions")
        else:
            return self.massfractions()
        
    @property
    def nelements(self):
        return len(self._elements)

    @property
    def ncompounds(self):
        return 1
        
    def markabsorber(self,symb=None,shells=None,fluolines=None,energybounds=None):
        """
        Args:
            symb(str): element symbol
        """
        for e in self._elements:
            if energybounds is None:
                energybounds2 = None
            else:
                energybounds2 = energybounds
            e.markabsorber(symb=symb,shells=shells,fluolines=fluolines,energybounds=energybounds2)

    def unmarkabsorber(self):
        for e in self._elements:
            e.unmarkabsorber()

    def hasabsorbers(self):
        return any([e.isabsorber for e in self._elements])

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
    
    def markinfo(self):
        yield "{}".format(self.name)
        yield " Scatterer: {}".format("yes" if self.isscatterer else "no")
        for e in self._elements:
            for s in e.markinfo():
                yield " {}".format(s)
        
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

        e_wfrac = self.massfractions()
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
                ret = sum(e_wfrac[e]*getattr(e,method)(E,environ=environ,**kwargs) for e in e_wfrac)

        return ret
        
    def refractive_index_delta(self,E,fine=False,decomposed=False,**kwargs):
        if hasattr(self,'structure') and fine:
            environ = self
        else:
            environ = None
        return self.refractive_index_delta_calc(E,self.massfractions(),self.density,environ=environ,**kwargs)
    
    def refractive_index_beta(self,E,fine=False,decomposed=False,**kwargs):
        if hasattr(self,'structure') and fine:
            environ = self
        else:
            environ = None
        return self.refractive_index_beta_calc(E,self.massfractions(),self.density,environ=environ,**kwargs)

    def get_energy(self,energyrange,defaultinc=1):
        """Get absolute energies (keV) from a relative energy range (eV)

        Args:
            energyrange(np.array): energy range relative to the absorber's edges (if any)
        """
        for e in self._elements:
            ret = e.get_energy(energyrange,defaultinc=defaultinc)
            if ret is not None:
                return ret
        return None

    def topymca(self,cfg,defaultthickness=1e-4):
        r = self.massfractions()
        massfractions = r.values()
        names = ['{}1'.format(e) for e in r]
        
        matname = self.pymcaname
        cfg["materials"][matname] = {'Comment': self.pymcacomment,
                                    'CompoundFraction': massfractions,
                                    'Thickness': defaultthickness,
                                    'Density': self.density,
                                    'CompoundList': names}
        return matname

    def tofisx(self,cfg,defaultthickness=1e-4):
        r = self.massfractions()
        massfractions = r.values()
        names = ['{}1'.format(e) for e in r]
        
        matname = self.pymcaname
        o = fisx.Material(matname, self.density, defaultthickness, self.pymcacomment)
        o.setCompositionFromLists(names,massfractions)
        cfg.addMaterial(o,errorOnReplace=False)
        
        return matname
        

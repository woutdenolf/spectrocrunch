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

from .simul import with_simulmetaclass

from . import noisepropagation

from ..materials.compoundfromformula import compound as compound
from ..materials.mixture import mixture
from ..materials.types import fractionType

import numpy as np

from . import emspectrum

from .. import ureg


class Scintillator(with_simulmetaclass()):
    """
    Class representing an area scintillator
    """

    def __init__(self, thickness=None, material=None, nvisperkeV=None, visspectrum=None):
        """
        Args:
            thickness(num): thickness in micron
            material(spectrocrunch.materials.compound|spectrocrunch.materials.mixture): scintillator compoisition
            nvisperkeV(num): number of VIS photons generated per keV
            visspectrum(visspectrum.discrete): VIS visspectrum
        """
        self.required(material,"material")
        self.required(thickness,"thickness")
        self.required(nvisperkeV,"nvisperkeV")
        self.required(visspectrum,"visspectrum")
        
        self.thickness = float(thickness)
        self.nvisperkeV = float(nvisperkeV)
        self.material = material
        self.visspectrum = visspectrum
        
    @staticmethod
    def doping(material,dopants,ftype):
        for el, frac in dopants.iteritems():
            material.addelement(el,frac,ftype)

    def absorption(self,energy):
        return 1-np.exp(-self.material.density*self.thickness*1e-4*self.material.mass_abs_coeff(energy))

    def attenuation(self,energy):
        return 1-self.transmission(energy)

    def transmission(self,energy):
        return np.exp(-self.material.density*self.thickness*1e-4*self.material.mass_att_coeff(energy))

    def propagate(self,N,energy):
        """Error propagation of a number of photons.
               
        Args:
            N(num|array): incomming number of photons with uncertainties
            energy(num|array): associated energies

        Returns:
            numpy.array,dict
        """
        
        # Absorption of X-rays
        probsuccess = self.absorption(energy)
        process = noisepropagation.bernouilli(probsuccess)
        Nout = noisepropagation.compound(N,process)

        # Fluorescence of visible photons
        # https://doi.org/10.1088/0031-9155/57/15/4885
        # http://arizona.openrepository.com/arizona/handle/10150/577317
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4669903/
        gain = energy*self.nvisperkeV
        process = noisepropagation.poisson(gain)
        Nout = noisepropagation.compound(Nout,process)
        
        return Nout
        
    def get_nrefrac(self):
        return self.material.nrefrac
    
class GGG_ID21(Scintillator):
    """
    Eu doped GGG
    """
    aliases = ["GGG ID21"]

    def __init__(self,thickness=None):
        """
        Args:
            thickness(num): thickness in micron
        """
        material = compound(["Gd","Ga","O"],[3,5,12],fractionType.mole,7.08,nrefrac=1.8,name="GGG")
        Scintillator.doping(material,{"Eu":0.03},fractionType.weight)
        
        visspectrum = emspectrum.discrete(ureg.Quantity([595,610,715],"nm"))
        
        super(GGG_ID21, self).__init__(thickness=thickness,material=material,nvisperkeV=32,visspectrum=visspectrum)

class LSO_ID21(Scintillator):
    """
    Tb doped LSO
    """
    aliases = ["LSO ID21"]

    def __init__(self,thickness=None):
        """
        Args:
            thickness(num): thickness in micron
        """
        material = compound(["Lu","Si","O"],[2,1,5],fractionType.mole,7.4,nrefrac=1.82,name="LSO")
        Scintillator.doping(material,{"Tb":0.03},fractionType.weight)
        
        visspectrum = emspectrum.discrete(ureg.Quantity(550,"nm"))
        
        super(LSO_ID21, self).__init__(thickness=thickness,material=material,nvisperkeV=40,visspectrum=visspectrum)

classes = Scintillator.clsregistry
aliases = Scintillator.aliasregistry
factory = Scintillator.factory



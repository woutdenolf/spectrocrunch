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

from ..common.classfactory import FactoryBase
import collections

from . import noisepropagation

from ..materials.compoundfromformula import compound as compound
from ..materials.mixture import mixture
from ..materials.types import fractionType

from .constants import LambdaEnergy

import numpy as np

class Scintillator(FactoryBase):
    """
    Class representing an area scintillator
    """
    registry = collections.OrderedDict()
    registry2 = collections.OrderedDict()

    def __init__(self, thickness=None, material=None, nvisperkeV=None):
        """
        Args:
            thickness(num): thickness in micron
            material(spectrocrunch.materials.compound|spectrocrunch.materials.mixture): scintillator compoisition
            nvisperkeV(num): number of VIS photons generated per keV
        """
        if thickness is None:
            raise RuntimeError("Thickness not defined for {}".format(self.__class__.__name__))
        if material is None:
            raise RuntimeError("Material not defined for {}".format(self.__class__.__name__))
        if nvisperkeV is None:
            raise RuntimeError("nvisperkeV not defined for {}".format(self.__class__.__name__))

        self.thickness = float(thickness)
        self.nvisperkeV = float(nvisperkeV)
        self.material = material

    @staticmethod
    def doping(material,dopants):
        for el, wfrac in dopants.iteritems():
            material.addelement(el,wfrac,fractionType.weight)

    def visenergyprofile(self,energy):
        """Energy profile of visible photons

        Args:
            energy(array-like): energy in keV

        Returns:
            np.array: normalized intensity
        """
        return energy*0

    def absorption(self,energy):
        return 1-np.exp(-self.material.density*self.thickness*1e-4*self.material.mass_abs_coeff(energy))

    def attenuation(self,energy):
        return 1-self.transmission(energy)

    def transmission(self,energy):
        return np.exp(-self.material.density*self.thickness*1e-4*self.material.mass_att_coeff(energy))

    def propagate(self,N,energy):
        """Error propagation of a number of photons.
               
        Args:
            N(num or numpy.array(uncertainties.core.Variable)): incomming number of photons with uncertainties
            energy(num or numpy.array): associated energies

        Returns:
            uncertainties.core.Variable or numpy.array(uncertainties.core.Variable)
        """

        # Absorption of X-rays
        probsuccess = self.absorption(energy)
        process = noisepropagation.bernouilli(probsuccess)
        Nout = noisepropagation.propagate(N,process)

        # Fluorescence of visible photons
        gain = energy*self.nvisperkeV
        process = noisepropagation.poisson(gain)
        Nout = noisepropagation.propagate(Nout,process)

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
        Scintillator.doping(material,{"Eu":0.03})
        super(GGG_ID21, self).__init__(thickness=thickness,material=material,nvisperkeV=32)

    def visenergyprofile(self,energy):
        """Energy profile of visible photons

        Args:
            energy(array-like): energy in keV

        Returns:
            np.array: normalized intensity
        """

        if hasattr(energy,"__iter__"):
            n = len(energy)
        else:
            n = 1

        ret = np.zeros(n,dtype=float)
        for l in [595,610,715]: # Eu
            ret[energy==LambdaEnergy(l)] = 1/3.
        return ret

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
        Scintillator.doping(material,{"Tb":0.03})
        super(LSO_ID21, self).__init__(thickness=thickness,material=material,nvisperkeV=40)
        
    def visenergyprofile(self,energy):
        """Energy profile of visible photons

        Args:
            energy(array-like): energy in keV

        Returns:
            np.array: normalized intensity
        """

        if hasattr(energy,"__iter__"):
            n = len(energy)
        else:
            n = 1

        ret = np.zeros(n,dtype=float)
        ret[energy==LambdaEnergy(550)] = 1 # Tb
        return ret

registry = Scintillator.registry
factory = Scintillator.factory


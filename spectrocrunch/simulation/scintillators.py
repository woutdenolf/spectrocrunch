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

from six import with_metaclass

from ..materials.compoundfromformula import compound as compound
from ..materials.mixture import mixture
from ..materials.types import fractionType

from .constants import LambdaEnergy

import numpy as np

class ScintillatorMeta(type):
    """
    Metaclass used to register all scintillator classes inheriting from Scintillator
    """
    def __init__(cls, name, bases, dct):
        cls.registry[name.lower()] = cls
        super(ScintillatorMeta, cls).__init__(name, bases, dct)

class Scintillator(with_metaclass(ScintillatorMeta, object)):
    """
    Class representing an area scintillator
    """
    registry = {}

    @classmethod
    def factory(cls, name, thickness, dopants=None):
        """
        Args:
            name(str): name of the scintillator

        Returns:
            Scintillator
        """
        name = name.lower()
        if name in cls.registry:
            return cls.registry[name](thickness, dopants=dopants)
        else:
            raise RuntimeError("Scintillator {} is not one of the registered scintillators: {}".format(name, cls.registry.keys()))

    def __init__(self, thickness=0, material=None, nvisperkeV=1):
        """
        Args:
            thickness(num): thickness in micron
            material(spectrocrunch.materials.compound|spectrocrunch.materials.mixture): scintillator compoisition
            nvisperkeV(num): number of VIS photons generated per keV
        """

        self.thickness = float(thickness)
        self.nvisperkeV = float(nvisperkeV)
        if material is None:
            self.material = compound([],[],fractionType.mole,0,name="vacuum")
        else:
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

    def attenuation(self,energy):
        return 1-self.attenuation(energy)

    def transmission(self,energy):
        return np.exp(-self.material.density*self.thickness*1e-4*self.material.mass_att_coeff(energy))

class GGG(Scintillator):
    """
    doped GGG
    """

    def __init__(self,thickness,dopants=None):
        """
        Args:
            thickness(num): thickness in micron
            dopants(Optional(dict)): elements:wfraction
        """
        material = compound(["Gd","Ga","O"],[3,5,12],fractionType.mole,7.1,nrefrac=1.8,name="GGG")
        if dopants is not None:
            Scintillator.doping(material,dopants)
        super(GGG, self).__init__(thickness=10,material=material,nvisperkeV=32)

    def visenergyprofile(self,energy):
        """Energy profile of visible photons

        Args:
            energy(array-like): energy in keV

        Returns:
            np.array: normalized intensity
        """
        energyvis = LambdaEnergy(610)
        
        if hasattr(energy,"__iter__"):
            n = len(energy)
        else:
            n = 1

        ret = np.zeros(n,dtype=float)
        ret[energy==energyvis] = 1
        return ret

class LSO(Scintillator):
    """
    doped GGG
    """

    def __init__(self,thickness,dopants=None):
        """
        Args:
            thickness(num): thickness in micron
            dopants(Optional(dict)): elements:wfraction
        """
        material = compound(["Lu","Si","O"],[2,1,5],fractionType.mole,7.4,nrefrac=1.82,name="LSO")
        if dopants is not None:
            Scintillator.doping(material,dopants)
        super(LSO, self).__init__(thickness=10,material=material,nvisperkeV=40)
        
    def visenergyprofile(self,energy):
        """Energy profile of visible photons

        Args:
            energy(array-like): energy in keV

        Returns:
            np.array: normalized intensity
        """
        energyvis = LambdaEnergy(540)
        
        if hasattr(energy,"__iter__"):
            n = len(energy)
        else:
            n = 1

        ret = np.zeros(n,dtype=float)
        ret[energy==energyvis] = 1
        return ret


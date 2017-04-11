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

from . import noisepropagation

from ..materials.compoundfromformula import compound as compound
from ..materials.mixture import mixture
from ..materials.types import fractionType

import numpy as np

class MaterialMeta(type):
    """
    Metaclass used to register all material classes inheriting from Material
    """
    def __init__(cls, name, bases, dct):
        cls.registry[name.lower().replace(" ","_")] = cls
        super(MaterialMeta, cls).__init__(name, bases, dct)

class Material(with_metaclass(MaterialMeta, object)):
    """
    Class representing an area material
    """
    registry = {}

    @classmethod
    def factory(cls, name, thickness, dopants=None):
        """
        Args:
            name(str): name of the material

        Returns:
            Material
        """
        name = name.lower().replace(" ","_")
        if name in cls.registry:
            return cls.registry[name](thickness, dopants=dopants)
        else:
            raise RuntimeError("Material {} is not one of the registered materials: {}".format(name, cls.registry.keys()))

    def __init__(self, thickness=0, material=None):
        """
        Args:
            thickness(num): thickness in micron
            material(spectrocrunch.materials.compound|spectrocrunch.materials.mixture): material composition
        """

        self.thickness = float(thickness)
        if material is None:
            self.material = compound([],[],fractionType.mole,0,name="vacuum")
        else:
            self.material = material

        super(Material, self).__init__(pmftype=pmftypes.bernouilli)

    def attenuation(self,energy):
        return 1-self.material.transmission(energy)

    def transmission(self,energy):
        return np.exp(-self.material.density*self.thickness*1e-4*self.material.mass_att_coeff(energy))

    def propagate(self,N,energy):
        """Error propagation of a number of photons.
               
        Args:
            N(uncertainties.unumpy.uarray): incomming number of photons with uncertainties
            energy(np.array): associated energies

        Returns:
            uncertainties.unumpy.uarray
        """

        # Transmission of X-rays
        probsuccess = self.transmission(energy)
        process = noisepropagation.bernouilli(probsuccess)
        Nout = noisepropagation.propagate(N,process)

        return Nout

#import xraylib
#for c in xraylib.GetCompoundDataNISTList():
#    c = MaterialMeta(c, (), {})





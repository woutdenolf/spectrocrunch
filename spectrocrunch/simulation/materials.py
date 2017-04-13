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
import xraylib

from ..common import classfactory

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
    def factory(cls, name):
        """
        Args:
            name(str): name of the material

        Returns:
            Material
        """
        name = name.lower().replace(" ","_")
        if name in cls.registry:
            return cls.registry[name]()
        else:
            raise RuntimeError("Material {} is not one of the registered materials: {}".format(name, cls.registry.keys()))

    def __init__(self, material, thickness):
        """
        Args:
            material(spectrocrunch.materials.compound): material composition
            thickness(num): thickness in micron
        """

        self.thickness = float(thickness)
        self.material = material

    def attenuation(self,energy):
        return 1-self.material.transmission(energy)

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

        # Transmission of X-rays
        probsuccess = self.transmission(energy)
        process = noisepropagation.bernouilli(probsuccess)
        Nout = noisepropagation.propagate(N,process)

        return Nout

class vacuum(Material):
    """
    Vacuum
    """

    def __init__(self,thickness):
        compound([],[],fractionType.mole,0,name="vacuum")
        super(ultralene, self).__init__(material,thickness)

class ultralene(Material):
    """
    Ultralene
    """

    def __init__(self,thickness):
        data = xraylib.GetCompoundDataNISTByName("Kapton Polyimide Film")
        material = compound(data["Elements"],data["massFractions"],fractionType.weight,data["density"],name=data["name"])
        super(ultralene, self).__init__(material,thickness)



try:
    
    for c in xraylib.GetCompoundDataNISTList():
        c = c.lower().replace(" ","_")
        #c = type(c, (Material,), {"__init__":})

    #print Material.registry
except ImportError:
    pass

factory = Material.factory


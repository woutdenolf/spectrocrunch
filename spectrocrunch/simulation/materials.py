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
import xraylib

import numpy as np

class Material(FactoryBase):
    """
    Class representing an area material
    """
    registry = collections.OrderedDict()
    registry2 = collections.OrderedDict()

    def __init__(self, material=None, thickness=None):
        """
        Args:
            material(list(spectrocrunch.materials.compound|mixture)): material composition
            thickness(num): thickness in micron
        """
        if thickness is None:
            raise RuntimeError("Thickness not defined for {}".format(self.__class__.__name__))
        if material is None:
            raise RuntimeError("Material not defined for {}".format(self.__class__.__name__))

        if hasattr(thickness,"__iter__"):
            self.thickness = np.asarray(thickness,dtype=float)
        else:
            self.thickness = np.asarray([thickness],dtype=float)

        if hasattr(material,"__iter__"):
            self.material = material
        else:
            self.material = [material]

        self.nlayers = len(self.material)

    def fluo_prob(self,energy,layer):
        self.material[layer].xrf_cross_section_decomposed(energy)
        
        return 1-self.transmission_prob(energy,layer)

    def transmission_prob(self,energy,layer):
        return np.exp(-self.material[layer].density*self.thickness[layer]*1e-4*self.material[layer].mass_att_coeff(energy))

    def propagate(self,N,energy):
        """Error propagation of a number of photons.
               
        Args:
            N(num or numpy.array(uncertainties.core.Variable)): incomming number of photons with uncertainties
            energy(num or numpy.array): associated energies

        Returns:
            uncertainties.core.Variable or numpy.array(uncertainties.core.Variable)
        """

        Nout = N

        if any([m.hasabsorbers() for m in self.material]):
            prob = self.fluo_prob
        else:
            prob = self.transmission_prob

        for layer in range(self.nlayers):
            probsuccess = prob(energy,layer)
            process = noisepropagation.bernouilli(probsuccess)
            Nout = noisepropagation.propagate(Nout,process)

        return Nout

class Vacuum(Material):
    """
    Vacuum
    """

    def __init__(self,thickness):
        material = compound([],[],fractionType.mole,0,name="vacuum")
        super(Vacuum, self).__init__(material=material,thickness=thickness)

class Ultralene(Material):
    """
    Ultralene
    """

    def __init__(self,thickness=None):
        data = xraylib.GetCompoundDataNISTByName("Kapton Polyimide Film")
        material = compound(data["Elements"],data["massFractions"],fractionType.weight,data["density"],name=data["name"])
        super(Ultralene, self).__init__(material=material,thickness=thickness)

def NistFactory(nistname):
    def __init__(self,thickness=None):
        data = xraylib.GetCompoundDataNISTByName(nistname)
        material = compound(data["Elements"],data["massFractions"],fractionType.weight,data["density"],name=data["name"])
        Material.__init__(self, material=material,thickness=thickness)
    newclass = type(nistname.replace(" ","_"), (Material,),{"__init__": __init__,"aliases":[nistname,nistname.lower()]})
    return newclass

try:
    for c in xraylib.GetCompoundDataNISTList():
        c = NistFactory(c)
except ImportError:
    pass

registry = Material.registry
factory = Material.factory



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

from .simul import SimulBase

from . import noisepropagation

from ..materials.compoundfromlist import compoundfromlist as compound
from ..materials.mixture import mixture
from ..materials.types import fractionType
from ..math.fit1d import lstsq

import numpy as np

class Multilayer(SimulBase):
    """
    Class representing an area material
    """
    
    def __init__(self, material=None, thickness=None):
        """
        Args:
            material(list(spectrocrunch.materials.compound|mixture)): material composition
            thickness(num): thickness in micron
        """
        self.required(material,"material")
        self.required(thickness,"thickness")
        
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

    def refinethickness(self,energy,absorbance,layerfixed=None):
        if layerfixed is None:
            layerfixed = []
        y = absorbance

        A = [self.material[layer].density*1e-4*self.material[layer].mass_att_coeff(energy) for layer in range(self.nlayers) if layer not in layerfixed]

        for layer in range(self.nlayers):
            if layer in layerfixed:
                y -= self.material[layer].density*1e-4*self.material[layer].mass_att_coeff(energy)

        A = np.vstack(A).T
        thickness = lstsq(A,y)
        ind = [layer not in layerfixed for layer in range(self.nlayers)]
        self.thickness[ind] = thickness

    def absorbance(self,energy):
        return sum([self.material[layer].density*self.thickness[layer]*1e-4*self.material[layer].mass_att_coeff(energy) for layer in range(self.nlayers)])

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
            prob = self.fluo_prob # Poisson instead of Bernouilli?
        else:
            prob = self.transmission_prob

        for layer in range(self.nlayers):
            probsuccess = prob(energy,layer)
            process = noisepropagation.bernouilli(probsuccess)
            Nout = noisepropagation.propagate(Nout,process)

        return Nout

classes = Multilayer.clsregistry
aliases = Multilayer.aliasregistry
factory = Multilayer.factory



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

import numpy as np

from ..patch.pint import ureg
from . import compound
from . import types
from ..sources import emspectrum
from ..utils import instance
from ..simulation.classfactory import with_metaclass
from ..math import noisepropagation


class Scintillator(with_metaclass()):
    """
    Class representing a scintillator
    """

    def __init__(self, thickness=None, material=None, nvisperkeV=None, visspectrum=None):
        """
        Args:
            thickness(num): thickness in micron
            material(Compound|Mixture): scintillator compoisition
            nvisperkeV(num): number of VIS photons generated per keV
            visspectrum(visspectrum.discrete): VIS visspectrum
        """
        self.thickness = thickness
        self.nvisperkeV = nvisperkeV
        self.material = material
        self.visspectrum = visspectrum

    def __getstate__(self):
        return {'thickness': self.thickness,
                'nvisperkeV': self.nvisperkeV,
                'material': self.material,
                'visspectrum': self.visspectrum}

    def __setstate__(self, state):
        self.thickness = state['thickness']
        self.nvisperkeV = state['nvisperkeV']
        self.material = state['material']
        self.visspectrum = state['visspectrum']

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.thickness == other.thickness and \
                self.nvisperkeV == other.nvisperkeV and \
                self.material == other.material and \
                self.visspectrum == other.visspectrum
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    @staticmethod
    def doping(material, dopants, ftype):
        material.addelements(list(dopants.keys()),
                             list(dopants.values()),
                             ftype)

    def absorption(self, energy):
        return 1-np.exp(-self.material.density*self.thickness*1e-4*self.material.mass_abs_coeff(energy))

    def attenuation(self, energy):
        return 1-self.transmission(energy)

    def transmission(self, energy):
        return np.exp(-self.material.density*self.thickness*1e-4*self.material.mass_att_coeff(energy))

    def propagate(self, N, energy, forward=True):
        """Error propagation of a number of photons.

        Args:
            N(num|array): incomming number of photons with uncertainties
            energy(num|array): associated energies

        Returns:
            numpy.array,dict
        """

        # Absorption of X-rays
        probsuccess = self.absorption(energy)

        # Fluorescence of visible photons
        # https://doi.org/10.1088/0031-9155/57/15/4885
        # http://arizona.openrepository.com/arizona/handle/10150/577317
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4669903/
        gain = energy*self.nvisperkeV

        N, probsuccess, gain = self.propagate_broadcast(N, probsuccess, gain)

        if instance.israndomvariable(N):
            if forward:
                proc1 = noisepropagation.bernouilli(probsuccess)
                proc2 = noisepropagation.poisson(gain)
            else:
                proc2 = noisepropagation.bernouilli(probsuccess)
                proc1 = noisepropagation.poisson(gain)

            Nout = noisepropagation.compound(N, proc1, forward=forward)
            Nout = noisepropagation.compound(Nout, proc2, forward=forward)
        else:
            if forward:
                Nout = N*(probsuccess*gain)
            else:
                Nout = N/(probsuccess*gain)

        return Nout

    def get_nrefrac(self):
        return self.material.nrefrac


class GGG_ID21(Scintillator):
    """
    Eu doped GGG
    """
    aliases = ["GGG ID21"]

    def __init__(self, thickness=None):
        """
        Args:
            thickness(num): thickness in micron
        """
        material = compound.Compound(["Gd", "Ga", "O"], [
                                     3, 5, 12], types.fraction.mole, 7.08, nrefrac=1.8, name="GGG")
        Scintillator.doping(material, {"Eu": 0.03}, types.fraction.mass)
        visspectrum = emspectrum.Discrete(
            ureg.Quantity([595, 610, 715], "nm"), [1, 1, 1])
        # http://www.esrf.eu/files/live/sites/www/files/Industry/documentation/F2_Scintillators.pdf
        super(GGG_ID21, self).__init__(thickness=thickness,
                                       material=material, nvisperkeV=32, visspectrum=visspectrum)


class LSO_ID21(Scintillator):
    """
    Tb doped LSO
    """
    aliases = ["LSO ID21"]

    def __init__(self, thickness=None):
        """
        Args:
            thickness(num): thickness in micron
        """
        material = compound.Compound(["Lu", "Si", "O"], [
                                     2, 1, 5], types.fraction.mole, 7.4, nrefrac=1.82, name="LSO")
        Scintillator.doping(material, {"Tb": 0.03}, types.fraction.mass)
        visspectrum = emspectrum.Discrete(ureg.Quantity(550, "nm"), 1)
        # http://www.esrf.eu/files/live/sites/www/files/Industry/documentation/F2_Scintillators.pdf
        super(LSO_ID21, self).__init__(thickness=thickness,
                                       material=material, nvisperkeV=40, visspectrum=visspectrum)


factory = Scintillator.factory
registry = Scintillator.clsregistry

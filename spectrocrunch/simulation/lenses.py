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

import numpy as np

class Lens(SimulBase):
    """
    Class representing a lens
    """

    def __init__(self, magnification=None, NA=None, transmission=None):
        """
        Args:
            magnification(num): magnification
            NA(num): numerical aperture
            transmission(num): transmission
        """
        self.required(magnification,"magnification")
        self.required(NA,"NA")
        self.required(transmission,"transmission")
        
        self.magnification = float(magnification)
        self.NA = float(NA)
        self.transmission = float(transmission)

    def propagate(self,N,energy,nrefrac=None):
        """Error propagation of a number of photons.
               
        Args:
            N(num or numpy.array(uncertainties.core.Variable)): incomming number of photons with uncertainties
            energy(num or numpy.array): associated energies
            nrefrac(num): refraction index of the scintillator

        Returns:
            uncertainties.core.Variable or numpy.array(uncertainties.core.Variable)
        """

        if nrefrac is None:
            ValueError("Refractive index of the scintillator not specified.")

        # Transmission of X-rays
        probsuccess = self.transmission*(self.NA*self.magnification/(self.magnification+1.)/(2.*nrefrac**2))**2
        process = noisepropagation.bernouilli(probsuccess)
        Nout = noisepropagation.propagate(N,process)

        # Repeat with number of energies
        if not hasattr(Nout,"__iter__") and hasattr(energy,"__iter__"):
            Nout = np.repeat(Nout,len(energy))

        return Nout

class mitutoyoid21_10x(Lens):
    """
    Mitutoyo M Plan Apo 20x 0.42 f = 200 mm
    """
    aliases = ["Mitutoyo ID21 10x"]

    def __init__(self):
        super(mitutoyoid21_10x, self).__init__(magnification=10, NA=0.42, transmission=0.95)

class mitutoyoid21_20x(Lens):
    """
    Mitutoyo M Plan Apo HR 10x 0.42 f = 200 mm
    """
    aliases = ["Mitutoyo ID21 20x"]

    def __init__(self):
        super(mitutoyoid21_20x, self).__init__(magnification=20, NA=0.42, transmission=0.95)

classes = Lens.clsregistry
aliases = Lens.aliasregistry
factory = Lens.factory


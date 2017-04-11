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

class LensMeta(type):
    """
    Metaclass used to register all lens classes inheriting from Lens
    """
    def __init__(cls, name, bases, dct):
        cls.registry[name.lower().replace(" ","_")] = cls
        super(LensMeta, cls).__init__(name, bases, dct)

class Lens(with_metaclass(LensMeta, object)):
    """
    Class representing an area lens
    """
    registry = {}

    @classmethod
    def factory(cls, name):
        """
        Args:
            name(str): name of the lens

        Returns:
            Lens
        """
        name = name.lower().replace(" ","_")
        if name in cls.registry:
            return cls.registry[name]()
        else:
            raise RuntimeError("Lens {} is not one of the registered lenses: {}".format(name, cls.registry.keys()))

    def __init__(self, magnification=1, NA=1, transmission=1):
        """
        Args:
            magnification(num): magnification
            NA(num): numerical aperture
            transmission(num): transmission
        """

        self.magnification = float(magnification)
        self.NA = float(NA)
        self.transmission = float(transmission)

    def propagate(self,N,energy,nrefrac):
        """Error propagation of a number of photons.
               
        Args:
            N(uncertainties.unumpy.uarray): incomming number of photons with uncertainties
            energy(np.array): associated energies
            nrefrac: refraction index of the scintillator

        Returns:
            uncertainties.unumpy.uarray
        """

        # Transmission of X-rays
        probsuccess = self.transmission*(self.NA*self.magnification/(self.magnification+1.)/(2.*nrefrac**2))**2
        process = noisepropagation.bernouilli(probsuccess)
        Nout = noisepropagation.propagate(N,process)

        return Nout

class mitutoyoid21_10x(Lens):
    """
    Mitutoyo M Plan Apo 20x 0.42 f = 200 mm
    """

    def __init__(self):
        super(mitutoyoid21_10x, self).__init__(magnification=10, NA=0.42, transmission=0.95)

class mitutoyoid21_20x(Lens):
    """
    Mitutoyo M Plan Apo HR 10x 0.42 f = 200 mm
    """

    def __init__(self):
        super(mitutoyoid21_20x, self).__init__(magnification=20, NA=0.42, transmission=0.95)


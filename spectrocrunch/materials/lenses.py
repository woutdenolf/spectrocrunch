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

from ..simulation.classfactory import with_metaclass
from ..simulation import noisepropagation
from ..common import instance
from .visirlib import Material

import numpy as np

class Lens(with_metaclass()):
    """
    Class representing a lens
    """

    def __init__(self, magnification=None, NA=None, thickness=None, material=None):
        """
        Args:
            magnification(num): magnification
            NA(num): numerical aperture
            thickness(num): in cm
            transmission(Material): transmission
        """
        self.required(magnification,"magnification")
        self.required(NA,"NA")
        self.required(thickness,"thickness")
        self.required(material,"material")
        
        self.magnification = float(magnification)
        self.NA = float(NA)
        self.thickness = float(thickness)
        self.material = material

    def transmission(self,visspectrum):
        linatt = np.asarray(self.material.linear_attenuation_coefficient(visspectrum.lines))
        ratios = visspectrum.ratios
        return np.sum(ratios * np.exp(-linatt*self.thickness))
    
    def lightyield(self,nrefrac):
        #air = visirlib.Material("other","air","Ciddor")
        #nmedium = np.mean(air.refractive_index(visspectrum.energies))
        nmedium = 1 # vacuum
        return ( np.tan(np.arcsin(self.NA/nmedium))*self.magnification/(2*(self.magnification+1.)*nrefrac) )**2
        
    def propagate(self,N,visspectrum,nrefrac=None,forward=True):
        """Error propagation of a number of photons.
               
        Args:
            N(unumpy.uarray): incomming number of photons with uncertainties
            visspectrum(emspectrum): visible light spectrum
            nrefrac(num): refraction index of the scintillator

        Returns:
            numpy.array
        """

        if nrefrac is None:
            ValueError("Refractive index of the scintillator not specified.")

        # Transmission of visible light
        #
        # http://onlinelibrary.wiley.com/doi/10.1118/1.598055/pdf
        # Non-Lambertian:
        # coupling efficiency = T.(M / (4.F#.(1+M).nscint))^2
        #
        #   1/So + 1/Si = 1/f
        #   F# = f/d
        #   M = Si/So
        #   NA = nmedium*sin(theta)
        #   tan(theta) = d/(2.f) = 1/(2.F#)
        #   2.F# = 1/tan(asin(NA/nmedium))
        #
        #   So = distance between scintillator and lens
        #   Si = distance between lens and ccd
        #   f = focal length
        #   M = geometrical maginification
        #   d = lens diameter
        #   nmedium = refractive index of air
        
        probsuccess = self.transmission(visspectrum)
        
        N,probsuccess = self.propagate_broadcast(N,probsuccess)
        
        lightyield = self.lightyield(nrefrac)
        
        if noisepropagation.israndomvariable(N):
            if forward:
                proc1 = noisepropagation.bernouilli(probsuccess)
                proc2 = noisepropagation.poisson(lightyield)
            else:
                proc2 = noisepropagation.bernouilli(probsuccess)
                proc1 = noisepropagation.poisson(lightyield)

            Nout = noisepropagation.compound(N,proc1,forward=forward)
            Nout = noisepropagation.compound(Nout,proc2,forward=forward)
        else:
            if forward:
                Nout = N*(probsuccess*lightyield)
            else:
                Nout = N/(probsuccess*lightyield)
                
        return Nout

class mitutoyoid21_10x(Lens):
    """
    Mitutoyo M Plan Apo 20x 0.42 f = 200 mm
    """
    aliases = ["Mitutoyo ID21 10x"]

    def __init__(self):
        super(mitutoyoid21_10x, self).__init__(magnification=10, NA=0.42, thickness=7.5, material=Material("glass","BK7","SCHOTT"))

class mitutoyoid21_20x(Lens):
    """
    Mitutoyo M Plan Apo HR 10x 0.42 f = 200 mm
    """
    aliases = ["Mitutoyo ID21 20x"]

    def __init__(self):
        super(mitutoyoid21_20x, self).__init__(magnification=20, NA=0.42, thickness=8., material=Material("glass","BK7","SCHOTT"))

classes = Lens.clsregistry
aliases = Lens.aliasregistry
factory = Lens.factory


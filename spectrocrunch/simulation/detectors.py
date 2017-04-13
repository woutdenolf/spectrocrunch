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

from uncertainties import ufloat
import numpy as np

class AreaDetectorMeta(type):
    """
    Metaclass used to register all detector classes inheriting from AreaDetector
    """
    def __init__(cls, name, bases, dct):
        cls.registry[name.lower().replace(" ","_")] = cls
        super(AreaDetectorMeta, cls).__init__(name, bases, dct)

class AreaDetector(with_metaclass(AreaDetectorMeta, object)):
    """
    Class representing an area detector
    """
    registry = {}

    @classmethod
    def factory(cls, name, etoadu, qe, aduoffset, darkcurrent, readoutnoise):
        """
        Args:
            name(str): name of the detector

        Returns:
            AreaDetector
        """
        name = name.lower().replace(" ","_")
        if name in cls.registry:
            return cls.registry[name](etoadu, qe, aduoffset, darkcurrent, readoutnoise)
        else:
            raise RuntimeError("Detector {} is not one of the registered detectors: {}".format(name, cls.registry.keys()))

    def __init__(self, etoadu, qe, aduoffset, darkcurrent, readoutnoise):
        """
        Args:
            etoadu(num): number of ADU per electron
            qe(num): detector quantum efficiency (e/ph)
            aduoffset(num): pixel intensity offset (ADU)
            darkcurrent(num): dark current (e/sec)
            readoutnoise(num): readout noise (e)
        """

        self.etoadu = float(etoadu)
        self.qe = float(qe)
        self.aduoffset = float(aduoffset)
        self.darkcurrent = ufloat(darkcurrent,np.sqrt(darkcurrent))
        self.readoutnoise = ufloat(0,readoutnoise)

    def propagate(self,N,energy,tframe=None,nframe=None):
        """Error propagation of a number of photons.
               
        Args:
            N(num or numpy.array(uncertainties.core.Variable)): incomming number of photons within a tframe timespan with uncertainties
            energy(num or numpy.array): associated energies
            tframe(num): time per frame (sec)
            nframe(num): number of frames (sec)

        Returns:
            uncertainties.core.Variable or numpy.array(uncertainties.core.Variable): detector signal in ADU
        """

        if tframe is None:
            ValueError("Frame exposure time not specified.")
        if nframe is None:
            ValueError("Number of frames not specified.")

        # Generation of electrons
        gain = self.qe
        process = noisepropagation.poisson(gain)
        Nout = noisepropagation.propagate(N,process) # units: e

        # Add dark current
        Nout += self.darkcurrent*tframe # units: e

        # Add read-out noise
        Nout += self.readoutnoise # units: e

        # Convert to ADU
        Nout *= self.etoadu # units: ADU

        # Add ADU offset
        Nout += self.aduoffset # units: ADU

        # Number of frames
        Nout *= nframe # units: ADU

        # Repeat with number of energies
        if not hasattr(Nout,"__iter__") and hasattr(energy,"__iter__"):
            Nout = np.repeat(Nout,len(energy))

        return Nout # units: ADU

class pcoedge55(AreaDetector):
    """
    PCO Edge 5.5
    """

    def __init__(self):
        super(pcoedge55, self).__init__(etoadu=1/0.45, qe=0.03, aduoffset=95.5, darkcurrent=7.4, readoutnoise=0.95)

factory = AreaDetector.factory


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

from uncertainties import ufloat

import numpy as np

class AreaDetector(SimulBase):
    """
    Class representing an area detector
    """

    def __init__(self, etoadu=None, qe=None, aduoffset=None, darkcurrent=None, readoutnoise=None):
        """
        Args:
            etoadu(num): number of ADU per electron
            qe(num): detector quantum efficiency (e/ph)
            aduoffset(num): pixel intensity offset (ADU)
            darkcurrent(num): dark current (e/sec)
            readoutnoise(num): readout noise (e)
        """
        self.required(etoadu,"etoadu")
        self.required(qe,"qe")
        self.required(aduoffset,"aduoffset")
        self.required(darkcurrent,"darkcurrent")
        self.required(readoutnoise,"readoutnoise")
        
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
            tframe(num or numpy.array): time per frame (sec)
            nframe(num or numpy.array): number of frames (sec)

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
        # https://spie.org/samples/PM170.pdf
        Nout += self.readoutnoise # units: e

        # Convert to ADU
        Nout *= self.etoadu # units: ADU

        # Add ADU offset
        Nout += self.aduoffset # units: ADU

        # Number of frames
        Nout = noisepropagation.repeat(nframe,Nout) # units: ADU

        # Repeat with number of energies
        if not hasattr(Nout,"__iter__") and hasattr(energy,"__iter__"):
            Nout = np.repeat(Nout,len(energy))

        return Nout # units: ADU

class pcoedge55(AreaDetector):
    """
    PCO Edge 5.5
    """
    aliases = ["PCO Edge 5.5"]

    def __init__(self):
        super(pcoedge55, self).__init__(etoadu=65536/30000., qe=0.7, aduoffset=95.5, darkcurrent=7.4, readoutnoise=0.95)

classes = AreaDetector.clsregistry
aliases = AreaDetector.aliasregistry
factory = AreaDetector.factory


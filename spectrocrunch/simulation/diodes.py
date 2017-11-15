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

from .simul import with_simulmetaclass

from . import noisepropagation

from ..detectors import diode


class Diode(with_simulmetaclass()):
    """
    Class representing a diode working on the following principle:

     Current generated from an incomming flux:
     I(A) = I(ph/s) . E(keV/ph)/Ehole(keV).1(e) . (1-exp(-mu.rho.d)) + D(A)

     Counts after an oscillator:
     I(cts/s) = Fmax(cts/s)/Vmax(V) . Rk(Ohm) . I(A) + F0(cts/s)

    """

    def __init__(self, pndiode=None):
        """
        Args:
            pndiode(PNdiode): charge generator
        """
        self.required(pndiode,"pndiode")
        
        self.pndiode = pndiode

    def setgain(self,Rout):
        """Set output resistance of the picoamperemeter(keithley)

        Args:
            Rout(num): output resistance (Ohm)

        Returns:
            num or array-like: C/ph
        """
        self.pndiode.setgain(Rout)

    def fluxtocurrent(self,energy,flux):
        """
        Args:
            flux(num or array-like): ph/s
            energy(num): keV

        Returns:
            num or array-like: current (A)
        """

        op = self.pndiode.op_fluxtocurrent(energy)
        return op(flux)

    def currenttoflux(self,energy,current):
        """
        Args:
            flux(num or array-like): A
            energy(num): keV

        Returns:
            num or array-like: flux (ph/s)
        """

        op = self.pndiode.op_currenttoflux(energy)
        return op(current)

    def cpstoflux(self,energy,cps):
        """
        Args:
            cps(num or array-like): cts/s
            energy(num): keV

        Returns:
            num or array-like: flux (ph/s)
        """

        op = self.pndiode.op_cpstoflux(energy)
        return op(cps)

    def fluxtocps(self,energy,flux):
        """
        Args:
            flux(num or array-like): ph/s
            energy(num): keV

        Returns:
            num or array-like: cps (cts/s)
        """

        op = self.pndiode.op_fluxtocps(energy)
        return op(flux)

    def propagate(self,N,energy,tframe=None,nframe=None):
        """Error propagation of a number of photons.
               
        Args:
            N(unumpy.uarray): incomming number of photons with uncertainties
            energy(numpy.array): associated energies
            tframe(num|numpy.array): time per frame (sec)
            nframe(num|numpy.array): number of frames (sec)

        Returns:
            uncertainties.core.Variable or numpy.array(uncertainties.core.Variable): detector signal in DU
        """

        if tframe is None:
            ValueError("Frame exposure time not specified.")
        if nframe is None:
            ValueError("Number of frames not specified.")

        # TODO: 
        

        return Nout # units: DU

class sxmidet(Diode):
    """
    SXM idet (cts in transmission)
    """
    aliases = ["sxmidet"]

    def __init__(self,model=True):
        super(sxmidet, self).__init__(pndiode=diode.sxmidet())

class xrdpico1(Diode):
    """
    SXM pico1 (current in transmission)
    """
    aliases = ["xrdpico1"]

    def __init__(self,model=True):
        super(xrdpico1, self).__init__(pndiode=diode.xrdpico1())

class sxmiodet1(Diode):
    """
    SXM iodet (cts before the KB)
    """
    aliases = ["sxmiodet1"]

    def __init__(self,model=True):
        super(sxmiodet1, self).__init__(pndiode=diode.sxmiodet1())

class sxmiodet2(Diode):
    """
    SXM iodet (cts before the KB)
    """
    aliases = ["sxmiodet2"]

    def __init__(self,model=True):
        super(sxmiodet2, self).__init__(pndiode=diode.sxmiodet2())

        


classes = Diode.clsregistry
aliases = Diode.aliasregistry
factory = Diode.factory


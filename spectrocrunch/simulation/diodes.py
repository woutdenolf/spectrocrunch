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

from uncertainties import ufloat

import numpy as np

from ..materials.compoundfromformula import compoundfromformula as compound

import silx.math.fit as fit

from ..resources import resource_filename

from scipy import interpolate

from ..math.linop import linop

from ..math.fit1d import linfit

import warnings

from .. import ureg

from . import constants

class _lineardevice(object):

    def __init__(self,gain,offset):
        """
        Args:
            gain(num): units/ph
            offset(num): units/s
        """
        self.gain = float(gain)
        self.offset = float(offset)

    def op_fluxtoups(self):
        return linop(self.gain,self.offset)

    def op_upstoflux(self):
        return (self.op_fluxtoups())**(-1)

    def calibrate(self,flux,ups):
        self.gain,self.offset = linfit(flux,ups)


class _oscillator(object):

    def __init__(self,Fmax,F0,Vmax):
        """
        Args:
            Fmax(num): maximal frequency (Hz)
            F0(num): frequency offset (Hz)
            Vmax(num): voltage corresponding to Fmax (V)
        """
        self.Fmax = float(Fmax)
        self.F0 = float(F0)
        self.Vmax = float(Vmax)
        self.pndiode = None

    def link(self,pndiode):
        self.pndiode = pndiode
        pndiode.oscillator = self
        
    def _countspercoulomb(self):
        """Return counts-per-coulomb

        Args:
            None

        Returns:
            num: cts/C
        """
        return self.Fmax/self.Vmax*self.pndiode.Rout

    def _offset(self):
        """Return frequency offset

        Args:
            None

        Returns:
            num: cts/s
        """
        return self.F0

    def op_currenttocps(self):
        """Operator to convert current to counts-per-second

        Args:
            None

        Returns:
            spectrocrunch.math.linop: slope (cts/C), intercept(cts/s)
        """
        return linop(self._countspercoulomb(),self._offset())

    def op_cpstocurrent(self):
        """Operator to convert counts-per-second to current

        Args:
            None

        Returns:
            spectrocrunch.math.linop: slope (C/cts), intercept(C/s)
        """
        return (self.op_currenttocps())**(-1)

    def op_cpstoflux(self,energy):
        """Operator to convert counts-per-second to flux

        Args:
            energy(num): keV

        Returns:
            spectrocrunch.math.linop: slope (ph/cts), intercept(ph/s)
        """
        return self.op_cpstocurrent()*self.pndiode.op_currenttoflux(energy)

    def op_fluxtocps(self,energy):
        """Operator to convert flux to counts-per-second 

        Args:
            energy(num): keV

        Returns:
            spectrocrunch.math.linop: slope (cts/ph), intercept(cts/s)
        """
        return self.pndiode.op_fluxtocurrent(energy)*self.op_currenttocps()


class _pndiode(object):

    def __init__(self, material, Rout, darkcurrent):
        """
        Args:
            material(compound|mixture): material composition
            Rout(num): output resistance (Ohm)
            darkcurrent(num): C/s
        """
        self.material = material
        self.setgain(Rout)
        #self.darkcurrent = noisepropagation.poisson(darkcurrent)
        self.darkcurrent = darkcurrent
        self.oscillator = None

    def link(self,oscillator):
        self.oscillator = oscillator
        oscillator.pndiode = self

    def _attenuation(self,energy,thickness):
        """Calculate attenuation: 1-exp(-mu.rho.thickness)

        Args:
            energy(num): keV
            thickness: micron

        Returns:
            num or array-like: attenuation
        """
        return (1-np.exp(-self.material.mass_att_coeff(energy)*self.material.density*(thickness*1e-4)))
    
    def _chargeperphoton(self,energy):
        """Return charge-per-photon generated for one or more photon energies

        Args:
            energy(num): keV

        Returns:
            num or array-like: C/ph
        """
        e = ureg.Quantity(1,ureg.elementary_charge).to(ureg.C).magnitude
        return e/self.ehole*energy*self._attenuation(energy,self.thickness)

    def setgain(self,Rout):
        """Set output resistance of the picoamperemeter(keithley)

        Args:
            Rout(num): output resistance (Ohm)

        Returns:
            None
        """
        self.Rout = float(Rout)

    def _offset(self):
        """Return dark current

        Args:
            None

        Returns:
            num: C/s
        """
        return self.darkcurrent

    def op_fluxtocurrent(self,energy):
        """Operator to convert flux to current

        Args:
            energy(num): keV

        Returns:
            spectrocrunch.math.linop: slope (C/ph), intercept(C/s)
        """
        return linop(self._chargeperphoton(energy),self._offset())

    def op_currenttoflux(self,energy):
        """Operator to convert current to flux

        Args:
            energy(num): keV

        Returns:
            spectrocrunch.math.linop: slope (ph/C), intercept(ph/s)
        """
        return (self.op_fluxtocurrent(energy))**(-1)
        
    def op_currenttocps(self):
        """Operator to convert current to counter-per-second

        Args:
            None

        Returns:
            spectrocrunch.math.linop: slope (cts/C), intercept(cts/s)
        """
        return self.oscillator.op_currenttocps()

    def op_cpstocurrent(self):
        """Operator to convert counter-per-second to current

        Args:
            None

        Returns:
            spectrocrunch.math.linop: slope (C/cts), intercept(C/s)
        """
        return self.oscillator.op_cpstocurrent()

    def op_cpstoflux(self,energy):
        """Operator to convert counts-per-second to flux

        Args:
            energy(num): keV

        Returns:
            spectrocrunch.math.linop: slope (ph/cts), intercept(ph/s)
        """
        return self.oscillator.op_cpstoflux(energy)

    def op_fluxtocps(self,energy):
        """Operator to convert flux to counts-per-second 

        Args:
            energy(num): keV

        Returns:
            spectrocrunch.math.linop: slope (cts/ph), intercept(cts/s)
        """
        return self.oscillator.op_fluxtocps(energy)


class _absolute_pndiode(_pndiode):

    def __init__(self, material, energy, response, model=True):
        """
        Args:
            material(list(spectrocrunch.materials.compound|mixture)): material composition
            energy(array-like): keV
            response(array-like): spectral responsivity (A/W)
            model(Optional(bool)): model the response for interpolation
        """
        super(_absolute_pndiode,self).__init__(material,1,0)

        self.model = model
        if self.model:
            self._fit(energy,response)
        else:
            self.finterpol = interpolate.interp1d(energy,response)

    def _fmodel(self,energy,thickness,ehole):   
        return self._attenuation(energy,thickness)/ehole

    def _fit(self,energy,response):
        """Calculate d and Ehole by fitting:

            I(A) = I(ph/s) . E(eV/ph) . e(C) . (1-exp(-mu.rho.d))/Ehole(eV) + D(A)
            P(W) = I(ph/s) . E(eV/ph) . e(J/eV)
            response(A/W) = (I(A)-D(A))/P(W) = (1-exp(-mu.rho.d))/Ehole(eV)
            response(mA/W) = (1-exp(-mu.rho.d)) / Ehole(keV)

        Args:
            energy(array-like): keV
            response(array-like): mA/W

        Returns:
            None
        """
        ehole = 1./np.max(response)
    
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            thickness = -np.log(1-response*ehole)/(self.material.mass_att_coeff(energy)*self.material.density)*1e4
            thickness = np.median(thickness[np.isfinite(thickness)])

        p, cov_matrix = fit.leastsq(self._fmodel, energy, response, [thickness,ehole])
        self.thickness,self.ehole = tuple(p)

        #import matplotlib.pyplot as plt
        #plt.figure()
        #plt.plot(energy,response,'o')
        #plt.plot(energy,self._fmodel(energy,thickness,ehole))
        #plt.plot(energy,self._fmodel(energy,self.thickness,self.ehole))
        #plt.show()

    def spectral_responsivity(self,energy):
        """Return spectral responsivity

        Args:
            energy(num or array-like): keV

        Returns:
            num or array-like: mA/W
        """
        if self.model:
            return self._fmodel(energy,self.thickness,self.ehole)
        else:
            return self.finterpol(energy)


class _calibrated_pndiode(_pndiode):

    def __init__(self, material, darkcurrent, energy, ratio, absdiode, Rout, model=True):
        """
        Args:
            material(list(spectrocrunch.materials.compound|mixture)): material composition
            darkcurrent(num): dark current (A)
            energy(num): keV
            ratio(array-like): diode current (dark subtracted) divided by the corresponding absdiode current at the same energy
            absdiode(_absolutediode): absolute diode used for calibration
            Rout(num): output resistance (Ohm)
            model(Optional(bool)): model the ratio for interpolation
        """
        super(_calibrated_pndiode,self).__init__(material,Rout,darkcurrent)

        self.model = model
        if self.model:
            self._fit(energy,ratio)
        else:
            self.finterpol = interpolate.interp1d(energy,ratio)

        self.absdiode = absdiode

    def _fmodel(self,energy,b,m1,m2):  
        return b+m1*np.exp(m2*energy)

    def _fit(self,energy,ratio):
        """
        Args:
            energy(array-like): keV
            ratio(array-like): (diode(A)-dark(A))/absdiode(A)

        Returns:
            None
        """
        # Elementary charge: e = 1.6e-19 C = 1.6e-19 J/eV
        # 1 eV = 1.6e-19 J
        # 1 A/W = 1 C/J = 1.6e-19 C/eV

        imin = np.argmin(energy)
        imax = np.argmax(energy)

        b = np.min(ratio)*0.95
        m2 = np.log((ratio[imax]-b)/(ratio[imin]-b))/(energy[imax]-energy[imin])
        m1 = (ratio[imax]-b)/np.exp(m2*energy[imax])

        p, cov_matrix = fit.leastsq(self._fmodel, energy, ratio, [b,m1,m2])
        self.b,self.m1,self.m2 = tuple(p)

        #import matplotlib.pyplot as plt
        #plt.figure()
        #plt.plot(energy,ratio,'o')
        #plt.plot(energy,self._fmodel(energy,b,m1,m2))
        #plt.plot(energy,self._fmodel(energy,self.b,self.m1,self.m2))
        #plt.show()

    def _chargeperphoton(self,energy):
        """Return charge-per-photon generated for one or more photon energies

        Args:
            energy(num): keV

        Returns:
            num or array-like: C/ph
        """
        if self.model:
            ratio = self._fmodel(energy,self.b,self.m1,self.m2)
        else:
            ratio = self.finterpol(energy)

        refresponse = self.absdiode.spectral_responsivity(energy)

        e = ureg.Quantity(1,ureg.elementary_charge).to(ureg.C).magnitude
        return e*energy*refresponse*ratio

class _noncalibrated_pndiode(_pndiode):

    def __init__(self, material, darkcurrent, thickness, ehole, Rout):
        """
        Args:
            material(list(spectrocrunch.materials.compound|mixture)): material composition
            darkcurrent(num): dark current (A)
            thickness(num): thickness in micron
            ehole(num): energy needed to create an electron-hole pair (keV)
            Rout(num): output resistance (Ohm)
        """
        super(_noncalibrated_pndiode,self).__init__(material,Rout,darkcurrent)

        self.thickness = float(thickness)
        self.ehole = float(ehole)


class Diode(with_simulmetaclass()):
    """
    Class representing a diode working on the following principle:

     Current generated from an incomming flux:
     I(A) = I(ph/s) . E(keV/ph)/Ehole(keV).e(C) . (1-exp(-mu.rho.d)) + D(A)

     Counts after an oscillator:
     I(cts/s) = F0(cts/s) + Fmax(cts/s)/Vmax(V) . Rk(Ohm) . I(A)

    """

    def __init__(self, pndiode=None):
        """
        Args:
            pndiode(_pndiode): charge generator
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

    def propagate(self,N,energy,tframe=None,nframe=None,withnoise=True):
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
        material = compound("Si",0,name="Si")

        ptb = np.loadtxt(resource_filename('id21/ptb.dat'))
        ptb = _absolute_pndiode(material, ptb[:,0], ptb[:,1],model=model)

        ird = np.loadtxt(resource_filename('id21/ird.dat'))
        ird = _calibrated_pndiode(material, 0, ird[:-1,0], ird[:-1,1], ptb, 1e5,model=model) # keithley: Rout=1e5 V/A (change with setgain)
        ird.link(_oscillator(1e6,0,10)) # VTOF: F0=0, Fmax=1e6, Vmax=10

        super(sxmidet, self).__init__(pndiode=ird)

class xrdpico1(Diode):
    """
    SXM pico1 (current in transmission)
    """
    aliases = ["xrdpico1"]

    def __init__(self,model=True):
        material = compound("Si",0,name="Si")

        pico1 = _noncalibrated_pndiode(material, 0, 1e3, constants.eholepair_si(), 10/2.1e-6, model=model) # darkcurrent = 0A, thickness = 1mm, ehole = ... keV, Rout = 10V / 2.1e-6A (change with setgain)

        super(xrdpico1, self).__init__(pndiode=pico1)

classes = Diode.clsregistry
aliases = Diode.aliasregistry
factory = Diode.factory


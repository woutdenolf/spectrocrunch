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

import numpy as np

from ..materials import compoundfromformula 

from ..materials import multilayer 

import silx.math.fit as fit

from ..resources import resource_filename

from scipy import interpolate

from ..math.linop import linop

from ..math.fit1d import linfit

import warnings

from .. import ureg

from ..common import units

from . import constants

#
#
# Unit of elementary charge: 1 e = 1.6e-19 C 
# Electronvolt: 1 eV = 1.6e-19 J
# 1 A/W = 1 C/J = 1.6e-19 C/eV = 1 e/eV
#
#
# Oscillator: current (A) to counts/sec (Hz)
#  I(Hz) = Fmax(Hz)/Vmax(V) . Vmax(V)/Imax(A) . I(A) + F0(Hz)
#    SXM: Vmax/Imax set
#    MICRODIFF: Imax set
#
# PNdiode: flux (ph/s) to current(A)
#  I(A) = I(ph/s).E(eV/ph).1(e).(1-T)/Ehole(eV) + D(A)
#       = P(W).1(e).(1-T)/Ehole(eV) + D(A)
#       = P(W).K(C/eV) + D(A)
#    T: transmission of the diode
#    P: radiative power
#    e: elementary charge
#
# An absolute diode measures both I(A) and P(W) so that
# another diode which only measures I(A) can be calibrated:
#  I(A) = P(W).FACT.K_abs(C/eV) + D(A)
#    KR = K(C/eV)/K_abs(C/eV)
#       = [I(A)-D(A)]/[I_abs(A)-D_abs(A)]
#       = (1-T)/Ehole(eV) . Ehole_abs(eV)/(1-T_abs)
#
# Upstream diode:
#  I0(ph/s) = I(ph/s).T.g
#    g: transmission of the optics (KB)
#
#

class Oscillator(object):
    """This class describes a voltage-to-frequency convertor
    """
    
    def __init__(self,Fmax=None,F0=None,Vmax=None):
        """
        Args:
            Fmax(num): maximal frequency (Hz)
            F0(num): frequency offset (Hz)
            Vmax(num): voltage corresponding to Fmax (V)
        """
        self.Fmax = Fmax.to("Hz")
        self.F0 = F0.to("Hz")
        self.Vmax = Vmax.to("volt")
        self.pndiode = None
        
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

    
class PNdiode(object):

    def __init__(self, material=None, Rout=None, darkcurrent=None, oscillator=None, secondarytarget=None):
        """
        Args:
            material(compound|mixture): material composition
            Rout(num): output resistance (Ohm)
            darkcurrent(num): C/s
            oscillator(Oscillator): 
            secondarytarget(multilayer): optional secondary target
        """
        self.material = material
        self.setgain(Rout)
        #self.darkcurrent = noisepropagation.poisson(darkcurrent)
        self.darkcurrent = darkcurrent
        self.secondarytarget = secondarytarget
        self.oscillator = oscillator
        if oscillator is not None:
            oscillator.pndiode = self
        
        # To be defined by the parent class for charge generation per photon
        self.thickness = None
        self.ehole = None

    def setgain(self,Rout):
        """Set output resistance of the picoamperemeter(keithley)

        Args:
            Rout(num): output resistance (Ohm)

        Returns:
            None
        """
        self.Rout = Rout.to("V/A")
        
    def _diode_attenuation(self,energy,thickness):
        """Calculate attenuation: 1-exp(-mu.rho.thickness)

        Args:
            energy(num): keV
            thickness: cm

        Returns:
            num or array-like: attenuation
        """
        return (1-np.exp(-self.material.mass_att_coeff(units.magnitude(energy,"keV"))*self.material.density*units.magnitude(thickness,"cm")))
    
    def _target_yield(self,energy):
        if self.secondarytarget is None:
            return energy,1
        else:
            return self.secondarytarget.spectrum()
        
    def _chargeperphoton(self,energy):
        """Return charge-per-photon generated for one or more photon energies

        Args:
            energy(num): keV

        Returns:
            num or array-like: C/ph
        """
        energy,targetyield = self._target_yield(energy)
        diodeatt = self._diode_attenuation(energy,self.thickness)
        chargeperphoton = (ureg.Quantity(1,ureg.elementary_charge)/self.ehole*units.Quantity(energy,"keV")).to("coulomb")
        return targetyield*diodeatt*chargeperphoton

    def _dark(self):
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
        return linop(self._chargeperphoton(energy),self._dark())

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
        return self.oscillator.op_cpstocurrent()*self.op_currenttoflux(energy)

    def op_fluxtocps(self,energy):
        """Operator to convert flux to counts-per-second 

        Args:
            energy(num): keV

        Returns:
            spectrocrunch.math.linop: slope (cts/ph), intercept(cts/s)
        """
        return self.op_fluxtocurrent(energy)*self.oscillator.op_currenttocps()


class AbsolutePNdiode(PNdiode):

    def __init__(self, energy=None, response=None, model=True, **kwargs):
        """
        Args:
            energy(array-like): keV
            response(array-like): spectral responsivity (A/W)
        """
        super(AbsolutePNdiode,self).__init__(**kwargs)

        response = response.to("A/W") # or e/eV

        self.model = model
        self._fit(energy,response)
        self.finterpol = interpolate.interp1d(energy,response)

    def spectral_responsivity(self,energy):
        """Return spectral responsivity

        Args:
            energy(num or array-like): keV

        Returns:
            num or array-like: A/W
        """
        if self.model:
            r = self._diode_attenuation(energy,self.thickness)/self.ehole*ureg.Quantity(1,ureg.elementary_charge)
        else:
            try:
                r = self.finterpol(energy)
            except:
                r = self._diode_attenuation(energy,self.thickness)/self.ehole*ureg.Quantity(1,ureg.elementary_charge)
         
        return units.Quantity(r,"A/W")

    def _fmodel(self,energy,thickness,ehole):
        return self._diode_attenuation(energy,thickness)/ehole

    def _fit(self,energy,response):
        """Calculate d and Ehole by fitting:

            I(A) = P(W) . 1(e) . (1-exp(-mu.rho.d))/Ehole(eV) + D(A)
            response(A/W) = (I(A)-D(A))/P(W) = 1(e).(1-exp(-mu.rho.d))/Ehole(eV)

        Args:
            energy(array-like): keV
            response(array-like): A/W

        Returns:
            None
        """

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            ehole = 1/np.max(response)
            muL = self.material.mass_att_coeff(energy.to("keV").magnitude)*self.material.density
            thickness = -np.log(1-(response*ehole).magnitude)/muL
            thickness = np.median(thickness[np.isfinite(thickness)])
            ehole = ehole.magnitude
        
        p, cov_matrix = fit.leastsq(self._fmodel, energy, response, [thickness,ehole])
        self.thickness = ureg.Quantity(p[0],"cm")
        self.ehole = ureg.Quantity(p[1],"eV")
        
        #import matplotlib.pyplot as plt
        #plt.figure()
        #plt.plot(energy,response,'o')
        #plt.plot(energy,self.spectral_responsivity(energy).to("A/W"))
        #ax = plt.gca()
        #ax.set_xlabel(energy.units)
        #ax.set_ylabel(response.units)
        #plt.show()


class CalibratedPNdiode(PNdiode):

    def __init__(self, energy=None, reponseratio=None, absdiode=None, model=True, **kwargs):
        """
        Args:
            energy(num): keV
            reponseratio(array-like): diode current (dark subtracted) divided by the corresponding absdiode current at the same energy
            absdiode(AbsolutePNdiode): absolute diode used for calibration
            model(Optional(bool)): model the ratio for interpolation
        """
        super(CalibratedPNdiode,self).__init__(**kwargs)

        self.model = model
        self._fit(energy,reponseratio)
        self.finterpol = interpolate.interp1d(energy,reponseratio)

        self.absdiode = absdiode

    def spectral_responsivity_ratio(self,energy):
        """Spectral responsivity ratio with the absolute diode

        Args:
            energy(num or array-like): keV

        Returns:
            num or array-like: 
        """
        if self.model:
            ratio = self._fmodel(energy,self.b,self.m1,self.m2)
        else:
            try:
                ratio = self.finterpol(energy)
            except:
                ratio = self._fmodel(energy,self.b,self.m1,self.m2)
        return ratio
        
    def _fmodel(self,energy,b,m1,m2):  
        return b+m1*np.exp(m2*units.magnitude(energy,"keV"))

    def _fit(self,energy,ratio):
        """
        Args:
            energy(array-like): keV
            ratio(array-like): (diode(A)-dark(A))/absdiode(A)

        Returns:
            None
        """

        energy = units.magnitude(energy,"keV")

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
        #plt.plot(energy,self.spectral_responsivity_ratio(energy))
        #plt.show()

    def _chargeperphoton(self,energy):
        """Return charge-per-photon generated for one or more photon energies

        Args:
            energy(num): keV

        Returns:
            num or array-like: C/ph
        """
        spectral_responsivity = self.absdiode.spectral_responsivity(energy)*self.spectral_responsivity_ratio(energy)
        return (units.Quantity(energy,"keV")*spectral_responsivity).to("coulomb")


class NonCalibratedPNdiode(PNdiode):

    def __init__(self, thickness=None, ehole=None, **kwargs):
        """
        Args:
            thickness(num): thickness in micron
            ehole(num): energy needed to create an electron-hole pair (keV)
        """
        super(NonCalibratedPNdiode,self).__init__(**kwargs)

        if thickness is None:
            self.thickness = None
        else:
            self.thickness = ureg.Quantity(thickness,"micrometer")
            
        if ehole is None:
            self.ehole = None
        else:
            self.ehole = ureg.Quantity(ehole,"eV")

    

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
        diodematerial = compoundfromformula.CompoundFromFormula("Si",0,name="Si")

        ptb = np.loadtxt(resource_filename('id21/ptb.dat'))
        energy = ureg.Quantity(ptb[:,0],"keV")
        response = ureg.Quantity(ptb[:,1],"milliampere/watt")
        ptb = AbsolutePNdiode(material=diodematerial,\
                            Rout=ureg.Quantity(1e5,"volt/ampere"),\
                            darkcurrent=ureg.Quantity(0,"ampere"),\
                            energy=energy,response=response,model=model)

        vtof = Oscillator(Fmax=ureg.Quantity(1e6,"Hz"),\
                        F0=ureg.Quantity(0,"Hz"),\
                        Vmax=ureg.Quantity(10,"volt"))

        ird = np.loadtxt(resource_filename('id21/ird.dat'))
        energy = ureg.Quantity(ird[:-1,0],"keV")
        responseratio = ird[:-1,1]
        ird = CalibratedPNdiode(material=diodematerial,\
                            Rout=ureg.Quantity(1e5,"volt/ampere"),\
                            darkcurrent=ureg.Quantity(0,"ampere"),\
                            energy=energy, reponseratio=responseratio, model=model ,absdiode=ptb,\
                            oscillator=vtof) 
        
        super(sxmidet, self).__init__(pndiode=ird)

class xrdpico1(Diode):
    """
    SXM pico1 (current in transmission)
    """
    aliases = ["xrdpico1"]

    def __init__(self,model=True):
        diodematerial = compoundfromformula.CompoundFromFormula("Si",0,name="Si")

        pico1 = NonCalibratedPNdiode(material=diodematerial,\
                                    Rout=ureg.Quantity(10,"volt")/ureg.Quantity(2.1e-6,"ampere"),\
                                    darkcurrent=ureg.Quantity(0,"ampere"),
                                    thickness=ureg.Quantity(1,"mm"),\
                                    ehole=constants.eholepair_si(),\
                                    model=model)

        super(xrdpico1, self).__init__(pndiode=pico1)

class sxmiodet1(Diode):
    """
    SXM iodet (cts before the KB)
    """
    aliases = ["sxmiodet1"]

    def __init__(self,model=True):

        diodematerial = compound("Si",0,name="Si")
        
        window = compoundfromname.compoundfromname("silicon nitride")
        coating = compoundfromformula.CompoundFromFormula("Ti",0,name="Ti")
        secondarytarget = multilayer.Multilayer(material=[coating,window],thickness=[1,0.5],anglein=0,angleout=135,solidangle=0.1*4*np.pi)

        vtof = Oscillator(Fmax=ureg.Quantity(1e6,"Hz"),\
                            F0=ureg.Quantity(0,"Hz"),\
                            Vmax=ureg.Quantity(10,"volt"))

        iodet1 = NonCalibratedPNdiode(material=diodematerial,\
                                Rout=ureg.Quantity(1e5,"volt/ampere"),\
                                darkcurrent=ureg.Quantity(0,"ampere"),\
                                thickness=ureg.Quantity(1,"mm"),\
                                ehole=constants.eholepair_si(),\
                                model=model,\
                                oscillator=vtof)
        
        super(sxmiodet1, self).__init__(pndiode=iodet1)

class sxmiodet2(Diode):
    """
    SXM iodet (cts before the KB)
    """
    aliases = ["sxmiodet2"]

    def __init__(self,model=True):
        diodematerial = compoundfromformula.CompoundFromFormula("Si",0,name="Si")
        
        window = compoundfromname.compoundfromname("silicon nitride")
        secondarytarget = multilayer.Multilayer(material=window,thickness=0.5,anglein=0,angleout=135,solidangle=0.1*4*np.pi)

        vtof = Oscillator(Fmax=ureg.Quantity(1e6,"Hz"),\
                        F0=ureg.Quantity(0,"Hz"),\
                        Vmax=ureg.Quantity(10,"volt"))

        iodet2 = NonCalibratedPNdiode(material=diodematerial,\
                        Rout=ureg.Quantity(1e5,"volt/ampere"),\
                        darkcurrent=ureg.Quantity(0,"ampere"),\
                        thickness=ureg.Quantity(1,"mm"),\
                        ehole=constants.eholepair_si(),\
                        model=model,\
                        oscillator = vtof,\
                        secondarytarget=secondarytarget)
        
        super(sxmiodet2, self).__init__(pndiode=iodet2)

        


classes = Diode.clsregistry
aliases = Diode.aliasregistry
factory = Diode.factory


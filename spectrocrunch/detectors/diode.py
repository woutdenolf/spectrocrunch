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
import warnings
import silx.math.fit as fit
import scipy.interpolate
import matplotlib.pyplot as plt
import logging
import scipy.optimize
import collections

from .. import ureg
from ..math import linop
from ..math import fit1d
from ..common import units
from ..materials import compoundfromformula 
from ..materials import compoundfromname
from ..materials import multilayer 
from ..materials import element
from ..resources import resource_filename
from ..optics import KB 
from ..common import constants
from ..math.common import round_sig
from ..geometries import diode as diodegeometries
from ..geometries import source as xraysources
from ..simulation.classfactory import with_metaclass
from ..simulation import noisepropagation
from ..common import instance
from . import base

logger = logging.getLogger(__name__)

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
#       = P(W).R(e/eV) + D(A)
#    T: transmission of the diode
#    P: radiative power
#    e: elementary charge
#    R: spectral responsivity
#
# Sample flux:
#  Is(ph/s) = I0(ph/s).Ts
#  I(ph/s) = I0(ph/s).Y = Is(ph/s).Y/Ts
#    I: flux seen by the diode
#    I0: the flux of the source
#    Is: flux seen by the sample
#    Y: yield of secondary target (includes solid angle of detector but not attenuation)
#    Ts: transmission of the source lines (optics,diode/target,filters)
#
# For multiple source lines (i) and multiple secondary lines (j):
#   Iis(ph/s) = I0(ph/s).wi.Tis
#   Is(ph/s) = sum_i[Iis] = I0(ph/s) . sum_i [wi.Tis]
#   Iij(ph/s) = I0(ph/s).Yij = Is(ph/s).Yij/sum_i[wi.Tis]
#
#   I(A) = SUM_ij [Iij(ph/s).Ej(eV/ph).1(e).(1-Tj)/Ehole(eV)] + D(A)
#        = SUM_ij [Iij(ph/s).Cj] + D(A)
#        = Is(ph/s) . SUM_i [SUM_j[Yij.Cj]] / SUM_k [wk.Tks] + D(A)
#        = Is(ph/s) . Cs + D(A)
#       
#   wi: source line fraction
#   Cj (C/ph): charge per photon hitting the diode
#   Cs (C/ph): charge per photon hitting the sample
#
# The source fraction at the sample position:
#   wis = Iis/sum_i[Iis] = wi.Tis / sum_i[wi.Tis]
#
# An absolute diode measures both I(A) and P(W) so that
# another diode which only measures I(A) can be calibrated:
#  I_abs(A) = P(W).R_abs(e/eV) + D_abs(A)
#  I(A) = P(W).FACT.R_abs(e/eV) + D(A)
#
#   FACT = R(e/eV)/R_abs(e/eV)
#        = [I(A)-D(A)]/[I_abs(A)-D_abs(A)] 
#        = (1-T)/Ehole(eV) . Ehole_abs(eV)/(1-T_abs)
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
    
    def tojson(self):
        kwargs =  {"Fmax":self.Fmax,"F0":self.F0,"Vmax":self.Vmax}
        return {"classname":self.__class__.__name__,"kwargs":kwargs}

    def __str__(self):
        return "y hertz = {:~e}/{:~} * x V + {:~}".format(self.Fmax,self.Vmax,self.F0)
    
    def _countspercoulomb(self):
        """Return counts-per-coulomb

        Args:
            None

        Returns:
            num: cts/C
        """
        return self.Fmax/self.Imax

    @property
    def Imax(self):
        return self.Vmax/self.pndiode.Rout

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
            callable: slope (cts/C), intercept(cts/s)
        """
        # These are the same:
        #return linop.Clip(cmin=None,cmax=self.Fmax) * linop.LinearOperator(self._countspercoulomb(),self._offset())
        #return linop.LinearOperator(self._countspercoulomb(),self._offset()) * linop.Clip(cmin=None,cmax=self.Imax)

        # Don't use clipping (xrfnormop may normalize iodet to a photon flux which may be outside the limit)
        return linop.LinearOperator(self._countspercoulomb(),self._offset())
        
    def op_cpstocurrent(self):
        """Operator to convert counts-per-second to current

        Args:
            None

        Returns:
            callable: slope (C/cts), intercept(C/s)
        """
        return self.op_currenttocps().inverse


class PNdiode(with_metaclass(base.SolidState)):

    ELCHARGE = ureg.Quantity(1,ureg.elementary_charge)
    #ELCHARGE = ureg.Quantity(1.6e-19,ureg.coulomb) # approx. in spec

    def __init__(self, Rout=None, darkcurrent=None, oscillator=None,\
                    secondarytarget=None, optics=None, beforesample=None,\
                    **kwargs):
        """
        Args:
            Rout(num): output resistance (Ohm)
            darkcurrent(num): C/s
            oscillator(Oscillator): 
            secondarytarget(multilayer): optional secondary target
        """
        self.setgain(Rout)
        self.setdark(darkcurrent)
        self.secondarytarget = secondarytarget
        self.beforesample = beforesample
        self.optics = optics
        self.oscillator = oscillator
        if oscillator is not None:
            oscillator.pndiode = self

        super(PNdiode,self).__init__(**kwargs)

    def generator(self):
        kwargs =  {"Rout":self.Rout.magnitude,"F0":self.F0,"Vmax":self.Vmax}
        return {"classname":self.__class__.__name__,"kwargs":kwargs}
        
    def __str__(self):
        fmt = "PN-diode:\n{}\n "\
              "Current:\n Gain = {:e}\n "\
              "Dark current = {}\n"\
              "Secondary target:\n {}\n"\
              "Optics:\n {}\n"\
              "Before sample: {}\n"\
              "Voltage-to-Frequency:\n {}"
                
        return fmt.format(super(PNdiode,self).__str__(),\
                self.Rout,self.darkcurrent,self.secondarytarget,\
                self.optics,self.beforesample,self.oscillator)    
    
    def setgain(self,Rout):
        """Set output resistance of the picoamperemeter(keithley)

        Args:
            Rout(num): output resistance (Ohm)

        Returns:
            None
        """
        self.Rout = units.Quantity(Rout,"V/A")
    
    def setdark(self,dark):
        """Set dark current

        Args:
            dark(num): dark current (Ampere)

        Returns:
            None
        """
        self.darkcurrent = units.Quantity(dark,"A")

    def _diode_absorbance(self,energy,thickness):
        # list self.absorbance but with customn thickness
        return self.material.mass_att_coeff(energy)*self.material.density*thickness
        
    def _diode_transmission(self,energy,thickness):
        # list self.transmission but with customn thickness
        return np.exp(-self._diode_absorbance(energy,thickness))

    def _diode_attenuation(self,energy,thickness):
        # list self.absorbance but with customn thickness
        return 1-self._diode_transmission(energy,thickness)
    
    def _spectral_responsivity_infthick(self):
        """Generated current per radiative power for an infinitely thick diode

        Args:
            energy(num): keV

        Returns:
            num or array-like: A/W or e/eV
        """
        return (self.ELCHARGE/self.ehole.to("eV")).to("A/W")
        
    def spectral_responsivity(self,energy):
        """Generated current per radiative power

        Args:
            energy(num): keV

        Returns:
            num or array-like: A/W or e/eV
        """
        return self.attenuation(energy)*self._spectral_responsivity_infthick()

    def plot_spectral_responsivity(self,energy):
        plt.plot(energy,self.spectral_responsivity(energy).to("A/W"))
        ax = plt.gca()
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Spectral responsivity (A/W)")
        ax.get_yaxis().get_major_formatter().set_useOffset(False) 
        plt.title(self.__class__.__name__)
        
        #ax.set_title("Current(A) = Responsivity(A/W) * Flux(ph/s) * Energy(eV) * e(J/eV)") # without secondary target, filters, ...
        #print self.fluxtocurrent(8.,1e8),self.spectral_responsivity(8.).magnitude*1e8*8e3*1.6e-19

    def plot_response(self,energy,flux,current=False):
        if current:
            y = [self.fluxtocurrent(en,flux).to("A").magnitude for en in energy]
            unit = "A"
        else:
            y = [self.fluxtocps(en,flux).to("Hz").magnitude for en in energy]
            unit = "Hz"
        
        plt.plot(energy,y)
        #plt.axhline(y=self.oscillator.Fmax.to("Hz").magnitude,label="max")
        ax = plt.gca()
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Response ({})".format(unit))
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        plt.title("{} @ {:.02e} ph/s, {:~.0e}".format(self.__class__.__name__,flux,self.Rout ))

    def _chargeperdiodephoton(self,energy):
        """Charge generated per photon absorbed by the diode: spectral responsivity multiplied by the energy

        Args:
            energy(num): keV

        Returns:
            num or array-like: C/ph
        """
        return (self.spectral_responsivity(energy)*units.Quantity(energy,"keV")).to("coulomb")
    
    def _transmission_optics(self,energy):
        """
        Args:
            energy(array): source energies in keV (dims: nE)
        
        Returns:
            
        """
        if self.optics is None:
            T = np.ones_like(energy)
        else:
            T = self.optics.transmission(energy)

        return T
        
    def _source_transmission(self,energy):
        """Transmission between sample and point-before-detection

        Args:
            energy(array): keV

        Returns:
            array: 
        """

        # Transmission of optics and filters-before-diode
        T = self._transmission_optics(energy)*\
            self.filter_transmission(energy,source=True)

        # Transmission of detection when before sample
        if self.beforesample:
            if self.secondarytarget is None:
                # diode itself in the beam
                T = T*self.transmission(energy)\
                     *self.filter_transmission(energy,source=False)
            else:
                # secondary target in the beam
                T = T*self.secondarytarget.transmission(energy)

        return T
    
    def _detection_rates(self,energy,weights):
        """Rate of line j due to line i (including source line weights and solid angle but not diode attenuation)
        
        Args:
            energy(array): source energies in keV (dims: nE)
            weights(array): source energies in keV (dims: nE)
            
        Returns:
            energy2(array): energies (keV) of lines detected by the diode (dims: nE2)
            rates(array): rates of lines detected by the diode (dims: nE x nE2)
            
        """
        
        if self.secondarytarget is None:
            # rates of the lines detected
            Y = weights[np.newaxis,:]
        else:
            nE = len(energy)
            
            # Spectrum generated from the target
            spectrum = self.secondarytarget.xrayspectrum(energy,weights=weights,withdetectorattenuation=False) 

            # Extract energies and rates (ph/phsource)
            kwargs = self.geometry.xrayspectrumkwargs()
            energy,Y = zip(*list(spectrum.spectrum(**kwargs)))
            
            energy = np.asarray(energy)
            nE2 = len(energy)
            Y = np.asarray(Y).reshape((nE,nE2))

        return energy,Y
    
    def _parse_source(self,energy,weights=None):
        energy = instance.asarray(energy)
        if weights is None:
            weights = np.ones_like(energy)
        else:
            weights = instance.asarray(weights)
        weights /= weights.sum() # make sure those are fractions
        return energy,weights
        
    def _chargepersamplephoton(self,energy,weights=None):
        """Charge generated per photon reaching the sample

        Args:
            energy(num|array): source energies in keV (dims: nE)
            weights(num|array): source line weights (dims: nE)

        Returns:
            num: C/ph
        """
        
        # Parse input
        energy,weights = self._parse_source(energy,weights=weights)
        
        # Transmission, rates and charge-per-diode-photon
        Ts = self._source_transmission(energy) # nE
        energy2,Yij = self._detection_rates(energy,weights) # nE2, nE x nE2
        Cj = self._chargeperdiodephoton(energy2) # nE2

        # I(A)  = Is(ph/s) . Cs(C) + D(A)
        # Cs(C) = SUM_i [SUM_j[Yij.Cj]] / SUM_k [wk.Tks]
        #
        # wi: source line fraction
        # Cj (C/ph): charge per photon hitting the diode
        # Cs (C/ph): charge per photon hitting the sample
        # Yij: rate of line j due to line i (including solid angle but not diode attenuation)

        # Sum over source lines: 
        return np.sum(Yij*Cj[np.newaxis,:])/np.sum(weights*Ts)
    
    def _calibrate_chargepersamplephoton(self,energy,Cs_calib,weights=None,caliboption=None,bound=False):
        """

        Args:
            energy(num|array): source energies in keV (dims: nE)
            Cs_calib(num): new charge-per-sample-photon (C)
            weights(num|array): source line weights (dims: nE)
            caliboption(str): "optics", "solidangle" or "thickness"
            
        Returns:
            None
        """

        # I(A) = Is(ph/s).Cs + D(A)
        #
        # Propagate cor to one of the component of Cs:
        #  Cs = SUM_i [wi.SUM_j[Yij.Cj]] / SUM_k [wk.Tks]
        #  Cj = (1-Tj).Cj_inf
        #  Cj_inf = Ej(eV/ph).1(e)/Ehole(eV)
        #  Yij = solidangle . ...
        #  Tj = 1-exp(-rho.mu(Ej).thickness)
        #
        #  Yij: rate of line j due to line i (including solid angle but not diode attenuation)
        #  Cj (C/ph): charge per photon hitting the diode
        #  Cs (C/ph): charge per photon hitting the sample

        if caliboption not in ["solidangle","thickness","optics"]:
            raise RuntimeError('caliboption must be "optics", "solidangle" or "thickness"')

        if caliboption=="solidangle":
            Cs = self._chargepersamplephoton(energy,weights=weights)
            
            # Yij' = Yij*Cs'/Cs
            sa = self.geometry.solidangle * units.magnitude(Cs_calib/Cs,"dimensionless")
            if sa<=0 or sa>(2*np.pi):
                logger.error("Diode solid angle of 4*pi*{} srad is not valid (possible wrong parameters: optics transmission, diode gain, diode thickness)".format(sa/(4*np.pi)))
            
            self.geometry.solidangle = sa
            Cscalc = self._chargepersamplephoton(energy,weights=weights)
 
        else:
            # TODO: some configurations have an analytical solution
            
            # Fluorescence does not change (expensive to calculate)
            energy,weights = self._parse_source(energy,weights=weights)
            energy2,Yij = self._detection_rates(energy,weights) # nE2, nE x nE2
            
            # Solve: Cs_calib-Cs(x) = 0
            if caliboption=="thickness":
                # x = diode thickness
                xorg = self.thickness
                if bound:
                    attenuationmax = 0.999
                    # use energy instead of energy2 because we need the maximal thickness
                    thicknessmax = -np.log(1-attenuationmax)/(self.material.mass_att_coeff(energy)*self.material.density)
                    x0 = (0.,np.max(thicknessmax))
                else:
                    x0 = xorg
                Cs = self._calibrate_thickness
                
            else:
                # x = transmission
                xorg = self.optics.transmission(energy)
                if len(energy)>1:
                    bound = False
                if bound:
                    x0 = (0.,1.)
                else:
                    x0 = xorg
                Cs = self._calibrate_optics

            # Solve: 1/Cs_calib - 1/Cs(x) = 0  (inverse because Cs are small values)
            func = lambda x: (1./Cs_calib-1./Cs(x,energy,weights,energy2,Yij)).to("1/coulomb").magnitude

            if bound:
                x,info = scipy.optimize.bisect(func,x0[0],x0[-1], full_output=True, disp=False)
                if not info.converged:
                    x = np.nan
            else:
                x,infodict,ier,emsg = scipy.optimize.fsolve(func, x0=x0, full_output=True)
                if ier!=1:
                    x = np.nan
                    
            # Check solution
            if np.isnan(x):
                logger.error("Diode calibration not succesfull")
                x = xorg

            if caliboption=="thickness":
                x = instance.asscalar(x)
                if x<=0:
                    logger.error("Diode thickness of {} um is not valid (possible wrong parameters: optics transmission, diode solid angle, diode gain)".format(x*1e-4))
                
                self.thickness = x
            else:
                if x<0 or x>1:
                    logger.error("Transmission of {} % is not valid (possible wrong parameters: diode thickness, solid angle or gain)".format(x*100))
                
                self.optics.set_transmission(energy,x)
            
            Cscalc = Cs(x,energy,weights,energy2,Yij)
   
        return Cscalc,Cs_calib
             
    def _calibrate_thickness(self,x,energy,weights,energy2,Yij):
        self.thickness = instance.asscalar(x)
        Ts = self._source_transmission(energy) # nE
        Cj = self._chargeperdiodephoton(energy2) # nE2
        return np.sum(Yij*Cj[np.newaxis,:])/np.sum(weights*Ts)
        
    def _calibrate_optics(self,x,energy,weights,energy2,Yij):
        self.optics.set_transmission(energy,x)
        Ts = self._source_transmission(energy) # nE
        Cj = self._chargeperdiodephoton(energy2) # nE2
        return np.sum(Yij*Cj[np.newaxis,:])/np.sum(weights*Ts)

    def samplelineweights(self,energy,weights):
        """Source weights after transmission
        
        Args:
            energy(num|array): source energies in keV (dims: nE)
            weights(num|array): source line weights (dims: nE)

        Returns:
            weights(num|array): source line weights at the sample position
        """
        
        # Parse input
        energy,weights = self._parse_source(energy,weights=weights)

        # wis = Iis/sum_i[Iis] = wi.Tis / sum_i[wi.Tis]
        weights *= self._source_transmission(energy)
        weights /= weights.sum()
        
        return weights
    
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
            callable: slope (C/ph), intercept(C/s)
        """
        return linop.LinearOperator(self._chargepersamplephoton(energy),self._dark())

    def op_currenttoflux(self,energy):
        """Operator to convert current to flux

        Args:
            energy(num): keV

        Returns:
            callable: slope (ph/C), intercept(ph/s)
        """
        return self.op_fluxtocurrent(energy).inverse
        
    def op_currenttocps(self):
        """Operator to convert current to counter-per-second

        Args:
            None

        Returns:
            callable: slope (cts/C), intercept(cts/s)
        """
        return self.oscillator.op_currenttocps()

    def op_cpstocurrent(self):
        """Operator to convert counter-per-second to current

        Args:
            None

        Returns:
            callable: slope (C/cts), intercept(C/s)
        """
        return self.oscillator.op_cpstocurrent()

    def op_cpstoflux(self,energy):
        """Operator to convert counts-per-second to flux

        Args:
            energy(num): keV

        Returns:
            callable: slope (ph/cts), intercept(ph/s)
        """
        return self.op_currenttoflux(energy)*self.op_cpstocurrent()

    def op_fluxtocps(self,energy):
        """Operator to convert flux to counts-per-second 

        Args:
            energy(num): keV

        Returns:
            callable: slope (cts/ph), intercept(cts/s)
        """
        return self.op_currenttocps()*self.op_fluxtocurrent(energy)
        
    def fluxtocurrent(self,energy,flux):
        """
        Args:
            flux(num or array-like): ph/s
            energy(num): keV

        Returns:
            num or array-like: current (A)
        """

        op = self.op_fluxtocurrent(energy)
        return op(units.Quantity(flux,"hertz")).to("ampere")

    def currenttoflux(self,energy,current):
        """
        Args:
            flux(num or array-like): A
            energy(num): keV

        Returns:
            num or array-like: flux (ph/s)
        """

        op = self.op_currenttoflux(energy)
        return op(units.Quantity(current,"ampere")).to("hertz")

    def cpstoflux(self,energy,cps):
        """
        Args:
            cps(num or array-like): cts/s
            energy(num): keV

        Returns:
            num or array-like: flux (ph/s)
        """
        op = self.op_cpstoflux(energy)
        return op(units.Quantity(cps,"hertz")).to("hertz")

    def fluxtocps(self,energy,flux):
        """
        Args:
            flux(num or array-like): ph/s
            energy(num): keV

        Returns:
            num or array-like: cps (cts/s)
        """
        op = self.op_fluxtocps(energy)
        return op(units.Quantity(flux,"hertz")).to("hertz")

    def currenttocps(self,current):
        """
        Args:
            current(num or array-like): A

        Returns:
            num or array-like: cps (cts/s)
        """
        op = self.op_currenttocps()
        return op(units.Quantity(current,"ampere")).to("hertz")
        
    def cpstocurrent(self,cps):
        """
        Args:
            cps(num or array-like): cts/s

        Returns:
            num or array-like: current (A)
        """
        op = self.op_cpstocurrent()
        return op(units.Quantity(cps,"hertz")).to("ampere")

    def xrfnormop(self,energy,expotime,reference,referencetime=None):
        """Operator to convert the raw diode signal to a flux normalizing signal.
           Usage: Inorm = I/op(iodet)
        
            XRF flux normalization:
            
            Measured XRF:
            Ixrf(cts) = F(cps).t(s).cxrf
            
            Desired XRF signal:
            Ixrf(cts)/norm = Fref(cps).tref(s).cxrf
            
            Normalization function to be apply on the raw diode signal Idiode
            norm = F.t/(Fref.tref) = cpstoflux(Idiode/t).t/(Fref.tref) = op(Idiode)
            op: x-> cpstoflux(x/t)/Fref.t/tref
            
            In case reference in counts instead of photons/sec
            Fref = round_sig(cpstoflux(Idioderef/t),2)
    
        Args:
            energy(num): keV
            expotime(num): sec
            reference(num): iodet (counts) or flux (photons/sec) to which the data should be normalized
            referencetime(Optional(num)): time to which the data should be normalized
            
        Returns:
            op(linop): raw diode conversion operator
            Fref(num): flux in photons/s to which the data is normalized after data/op(diode)
            tref(num): time in s to which the data is normalized after data/op(diode)
        """

        # Convert from counts to photons/sec
        # op: x-> cpstoflux(x/t)
        op = self.op_cpstoflux(energy)
        t = units.Quantity(expotime,"s")
        op.m /= t
        
        # Reference flux to which the XRF signal should be normalized
        if reference.units==ureg.hertz: # photons/sec
            Fref = reference
        elif reference.units==ureg.dimensionless: # counts
            # Fref = op(Idioderef)
            Fref = units.Quantity(round_sig(units.magnitude(op(reference),"hertz"),2),"hertz")
        else:
            raise RuntimeError("Reference {} should be in photons/sec (flux) or counts (iodet).".format(reference))
            
        # Convert from counts to counts at reference flux Fref
        op.m /= Fref
        op.b /= Fref
        if referencetime is not None:
            tref = units.Quantity(tref,"s")
            op.m *= t/tref
            op.b *= t/tref
        else:
            tref = t
            
        op.m = units.magnitude(op.m,"dimensionless")
        op.b = units.magnitude(op.b,"dimensionless")

        return op,Fref.to("hertz").magnitude,units.magnitude(tref,"s")
    
    def fluxop(self,energy,expotime):
        """Operator to convert the raw diode signal to a flux.
        
        Args:
            energy(num): keV
            expotime(num): sec
            
        Returns:
            op(linop): raw diode conversion operator
        """
        
        # Convert from counts to photons/sec
        # op: x-> cpstoflux(x/t)
        op = self.op_cpstoflux(energy)
        t = units.Quantity(expotime,"s")
        op.m /= t

        op.m = units.magnitude(op.m,"Hertz")
        op.b = units.magnitude(op.b,"Hertz")

        return op

    def gainfromcps(self,energy,cps,fluxest):
        """Estimate the gain, assuming it is 10^x V/A
        
        Args:
            energy(num): keV
            cps(num): 
            fluxest(num): estimated flux 

        Returns:
            None
        """
        
        with np.errstate(divide='ignore', invalid='ignore'):
            fluxcalc = self.cpstoflux(energy,cps)
            r = np.log10(units.magnitude(fluxcalc/units.Quantity(fluxest,"hertz"),"dimensionless"))
            r = int(round(np.nanmedian(r)))
            R = units.Quantity(10**(np.log10(units.magnitude(self.Rout,"V/A"))+r),"V/A")
            self.setgain(R)

    def darkfromcps(self,darkcps):
        """Dark current from cps
        
        Args:
            darkcps(num): measured dark counts per second

        Returns:
            None
        """
        darkcurrent = self.cpstocurrent(np.nanmedian(darkcps))
        self.setdark(darkcurrent)
    
    def fluxcpsinfo(self,energy):
        ret = collections.OrderedDict()
        Cs = self._chargepersamplephoton(energy).to("e")
        ret["Energy"] = "{} keV".format(energy)
        ret["Slope"] = "{:f} e$^-$/ph".format(Cs.to("e").magnitude)
        ret["Intercept"] = "{:e} e$^-$/s".format(self.darkcurrent.to("e/s").magnitude)
        return ret

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
        self.finterpol = scipy.interpolate.interp1d(energy,response,bounds_error=False,fill_value="extrapolate")

    def spectral_responsivity(self,energy):
        """Return spectral responsivity

        Args:
            energy(num or array-like): keV

        Returns:
            num or array-like: A/W
        """
        if self.model:
            r = self.attenuation(energy)*(self.ELCHARGE/self.ehole.to("eV"))
        else:
            try:
                r = self.finterpol(energy)
            except:
                r = self.attenuation(energy)*(self.ELCHARGE/self.ehole.to("eV"))
         
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
            muL = self.material.mass_att_coeff(energy)*self.material.density
            thickness = -np.log(1-(response*ehole).magnitude)/muL
            thickness = np.median(thickness[np.isfinite(thickness)])
            ehole = ehole.magnitude
        
        p, cov_matrix = fit.leastsq(self._fmodel, energy, response, [thickness,ehole])
        self.thickness = p[0] # "cm"
        self.ehole = units.Quantity(p[1],"eV") # eV


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
        fill_value = (reponseratio[0],reponseratio[-1])
        self.finterpol = scipy.interpolate.interp1d(energy,reponseratio,bounds_error=False,fill_value=fill_value)

        self.absdiode = absdiode

    def _spectral_responsivity_ratio(self,energy):
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
    
    def spectral_responsivity(self,energy):
        return self.absdiode.spectral_responsivity(energy)*self._spectral_responsivity_ratio(energy)
        
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

        imin = np.argmin(energy)
        imax = np.argmax(energy)

        b = np.min(ratio)*0.95
        m2 = np.log((ratio[imax]-b)/(ratio[imin]-b))/(energy[imax]-energy[imin])
        m1 = (ratio[imax]-b)/np.exp(m2*energy[imax])

        p, cov_matrix = fit.leastsq(self._fmodel, energy, ratio, [b,m1,m2])
        self.b,self.m1,self.m2 = tuple(p)
        
        
class NonCalibratedPNdiode(PNdiode):

    def __init__(self, **kwargs):
        super(NonCalibratedPNdiode,self).__init__(**kwargs)

    def calibrate(self,cps,sampleflux,energy,weights=None,caliboption=None,fixdark=False,fluxmin=0,fluxmax=np.inf):
        """Calibrate with another diode measuring the flux at the sample position

        Args:
            cps(array): count rate measured by this diode (cts/s)
            sampleflux(array): flux measured at the sample position (ph/s)
            energy(num|array): source lines (keV)
            weights(Optional(num|array)): source line weights
            caliboption(str): "optics", "solidangle" or "thickness"
            fixdark(Optional(num)): fix dark current
            
        Returns:
            None
        """

        # I(A) = Is(ph/s).slope + intercept
        #      = Is(ph/s).Cs + D(A)
        # Cs: charge per sample photon
        current = self.cpstocurrent(units.Quantity(cps,"hertz"))

        x = units.magnitude(sampleflux,"hertz")
        indfit = (x>=units.magnitude(fluxmin,"hertz")) & (x<=units.magnitude(fluxmax,"hertz"))
        
        if sum(indfit)<=2:
            raise RuntimeError("Not enough data points with a flux in [{:e},{:e}] (data range [{:e},{:e}])".format(fluxmin,fluxmax,np.min(x),np.max(x)))

        x = x[indfit]
        if fixdark:
            # Fit dark-subtracted current vs flux
            intercept = units.magnitude(self.darkcurrent,"ampere")
            y = units.magnitude(current,"ampere")
            y = y[indfit]
            slope = fit1d.linfit_zerointercept(x,y-intercept)
        else:
            # Fit current vs flux
            y = units.magnitude(current,"ampere")
            y = y[indfit]
            slope,intercept = fit1d.linfit(x,y)
            
            # Set dark current:
            self.setdark(intercept)
            
        # Correlation coefficient
        ycalc = intercept + slope*x
        R2 = 1-sum((y-ycalc)**2)/sum((y-np.mean(y))**2)
        
        # Set diode thickness, solid angle or transmission
        slope = units.Quantity(slope,"ampere/hertz")
        Cscalc,Cs_calib = self._calibrate_chargepersamplephoton(energy,slope,weights=weights,caliboption=caliboption)

        info = "Diode calibration:\n Energy: {} keV"\
               "\n Electron-hole pairs per photon hitting the sample: {:~} (experiment: {:~})"\
               "\n Electron-hole pairs per second (dark): {:~e} "\
               "\n R^2 = {}".\
                format(energy,Cscalc.to("e"),Cs_calib.to("e"),self.darkcurrent.to("e/s"),R2)
        logger.info(info)
        
        ret = self.fluxcpsinfo(energy)
        ret["R$^2$"] = R2
        return ret


class SXM_PTB(AbsolutePNdiode):

    aliases = ["ptb"]
    
    def __init__(self,**kwargs):
        attenuators = {}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":30e-4}
        kwargs["attenuators"] = attenuators
        
        ptb = np.loadtxt(resource_filename('id21/ptb.dat'))
        energy = ptb[:,0] # keV
        response = ureg.Quantity(ptb[:,1],"milliampere/watt")

        super(SXM_PTB,self).__init__(\
                        Rout=ureg.Quantity(1e5,"volt/ampere"),\
                        darkcurrent=ureg.Quantity(0,"ampere"),\
                        energy=energy,response=response,\
                        beforesample=False,**kwargs)

class SXM_IDET(CalibratedPNdiode):
    # Centronic OSD 50-3T
    
    aliases = ["idet"]
    
    def __init__(self,**kwargs):
    
        kwargs["attenuators"] = {}
        kwargs["attenuators"]["Detector"] = {"material":element.Element('Si'),"thickness":300e-4} # thickness not used because of PTB calibration
        kwargs["ehole"] = constants.eholepair_si() # ehole not used because of PTB calibration
        
        absdiode = SXM_PTB(model=True)

        vtof = Oscillator(Fmax=ureg.Quantity(1e6,"Hz"),\
                        F0=ureg.Quantity(32,"Hz"),\
                        Vmax=ureg.Quantity(10,"volt"))

        ird = np.loadtxt(resource_filename('id21/ird.dat'))
        energy = ird[:,0] # keV
        responseratio = ird[:,1]
        
        #energy = ird[:-1,0] # keV
        #responseratio = ird[:-1,1]

        super(SXM_IDET,self).__init__(\
                        Rout=ureg.Quantity(1e5,"volt/ampere"),\
                        darkcurrent=ureg.Quantity(1.08e-10,"ampere"),\
                        energy=energy,reponseratio=responseratio,\
                        absdiode=absdiode,oscillator=vtof,beforesample=False,**kwargs)
                        

class SXM_IODET1(NonCalibratedPNdiode):
    # International Radiation Detectors (IRD), AXUV-PS1-S
    # Keithley K428 (10V max analog output)
    # NOVA N101VTF voltage-to-frequency converter (Fmax=1e6, F0=0Hz)
    # P201 counter board

    aliases = ["iodet1"]

    def __init__(self,**kwargs):

        kwargs2 = {}
        kwargs2["anglein"] = kwargs.pop("anglein",62)
        kwargs2["angleout"] = kwargs.pop("angleout",49)
        kwargs2["azimuth"] = kwargs.pop("azimuth",0)
        kwargs2["solidangle"] = kwargs.pop("solidangle",4*np.pi*0.4)
        if "source" in kwargs:
            kwargs2["source"] = kwargs.pop("source")
        else:
            kwargs2["source"] = xraysources.factory("synchrotron")
        kwargs2["detector"] = self
        geometry = diodegeometries.factory("DiodeGeometry",**kwargs2)

        kwargs["attenuators"] = {}
        kwargs["attenuators"]["Detector"] = {"material":element.Element('Si'),"thickness":0.1}
        kwargs["ehole"] = constants.eholepair_si()
        
        window = compoundfromname.compoundfromname("silicon nitride")
        coating = element.Element('Ti')
        secondarytarget = multilayer.Multilayer(material=[coating,window],thickness=[500e-7,500e-7],geometry=geometry)
        
        optics = kwargs.pop("optics",True)
        if optics:
            optics = KB.KB()
        else:
            optics = None
        
        vtof = Oscillator(Fmax=ureg.Quantity(1e6,"Hz"),\
                            F0=ureg.Quantity(247,"Hz"),\
                            Vmax=ureg.Quantity(10,"volt"))
                        
        super(SXM_IODET1,self).__init__(\
                            Rout=ureg.Quantity(1e5,"volt/ampere"),\
                            darkcurrent=ureg.Quantity(0,"ampere"),\
                            oscillator=vtof,\
                            secondarytarget=secondarytarget,\
                            optics=optics,\
                            beforesample=True,\
                            **kwargs)
                        
class SXM_IODET2(NonCalibratedPNdiode):
    # International Radiation Detectors (IRD), AXUV-PS1-S
    # Keithley K428 (10V max analog output)
    # NOVA N101VTF voltage-to-frequency converter (Fmax=1e6, Vmax=0Hz)
    # P201 counter board
    
    aliases = ["iodet2"]

    def __init__(self,**kwargs):
    
        kwargs2 = {}
        kwargs2["anglein"] = kwargs.pop("anglein",62)
        kwargs2["angleout"] = kwargs.pop("angleout",49)
        kwargs2["azimuth"] = kwargs.pop("azimuth",0)
        kwargs2["solidangle"] = kwargs.pop("solidangle",4*np.pi*0.4)
        if "source" in kwargs:
            kwargs2["source"] = kwargs.pop("source")
        else:
            kwargs2["source"] = xraysources.factory("synchrotron")
        kwargs2["detector"] = self
        geometry = diodegeometries.factory("DiodeGeometry",**kwargs2)
        
        kwargs["attenuators"] = {}
        kwargs["attenuators"]["Detector"] = {"material":element.Element('Si'),"thickness":0.1}
        kwargs["attenuators"] = attenuators
        kwargs["ehole"] = constants.eholepair_si()
        
        window = compoundfromname.compoundfromname("silicon nitride")
        secondarytarget = multilayer.Multilayer(material=[window],thickness=[500e-7],geometry=geometry)
        
        optics = kwargs.pop("optics",True)
        if optics:
            optics = KB.KB()
        else:
            optics = None
            
        vtof = Oscillator(Fmax=ureg.Quantity(1e6,"Hz"),\
                            F0=ureg.Quantity(247,"Hz"),\
                            Vmax=ureg.Quantity(10,"volt"))
                        
        super(SXM_IODET2,self).__init__(\
                            Rout=ureg.Quantity(1e5,"volt/ampere"),\
                            darkcurrent=ureg.Quantity(0,"ampere"),\
                            oscillator=vtof,\
                            secondarytarget=secondarytarget,\
                            optics=optics,\
                            beforesample=True,\
                            **kwargs)

class XRD_IDET(NonCalibratedPNdiode):
    aliases = ["pico1"]
    
    def __init__(self,**kwargs):
    
        kwargs["attenuators"] = {}
        kwargs["attenuators"]["Detector"] = {"material":element.Element('Si'),"thickness":0.1}
        kwargs["ehole"] = constants.eholepair_si()

        super(XRD_IDET,self).__init__(\
                                Rout=ureg.Quantity(10,"volt")/ureg.Quantity(2.1e-6,"ampere"),\
                                darkcurrent=ureg.Quantity(0,"ampere"),beforesample=False,\
                                **kwargs)

factory = PNdiode.factory


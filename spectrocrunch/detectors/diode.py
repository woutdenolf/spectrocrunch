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
from scipy import interpolate

from .. import ureg
from ..math.linop import linop
from ..math.fit1d import linfit
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
#  Is(ph/s) = I0(ph/s).Ts.To
#  I(ph/s) = I0(ph/s).Y = Is(ph/s).Y/(Ts.To)
#    I: flux seen by the diode
#    I0: the flux of the source
#    Is: flux seen by the sample
#    Y: yield of secondary target (includes solid angle of detector)
#    Ts: transmission of the secondary target
#    To: transmission of the optics (e.g. reflectivity of the KB)
#
# For multiple source lines (i) and multiple secondary lines (j):
#   Iis(ph/s) = I0(ph/s).Fi.Tis.Tio
#   Is(ph/s) = sum_i[Iis] = I0(ph/s) . sum_i [Fi.Tis.Tio]
#   Iij(ph/s) = I0(ph/s).Fi.Yij = Is(ph/s).Fi/sum_i[Fi.Tis.Tio].Yij
#
#   I(A) = SUM_ij [Iij(ph/s).Ej(eV/ph).1(e).(1-Tj)/Ehole(eV)] + D(A)
#        = SUM_ij [Iij(ph/s).Cj] + D(A)
#        = Is(ph/s) . SUM_i [Fi.SUM_j[Yij.Cj]] / SUM_k [Fk.Tks.Tko] + D(A)
#        = Is(ph/s) . Cs + D(A)
#       
#   Fi: source line fraction
#   Cj (C/ph): charge per photon hitting the diode
#   Cs (C/ph): charge per photon hitting the sample
#
# The source fraction at the sample position:
#   Fis = Iis/sum_i[Iis] = Fi.Tis.Tio / sum_i[Fi.Tis.Tio]
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
    
    def __str__(self):
        return "y hertz = {}/{} * x volt + {}".format(self.Fmax,self.Vmax,self.F0)
    
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


class PNdiode(with_metaclass(base.SolidAngle)):

    ELCHARGE = ureg.Quantity(1,ureg.elementary_charge)
    #ELCHARGE = ureg.Quantity(1.6e-19,ureg.coulomb) # approx. in spec

    def __init__(self, Rout=None, darkcurrent=None, oscillator=None, secondarytarget=None, optics=None, **kwargs):
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
        self.optics = optics
        self.oscillator = oscillator
        if oscillator is not None:
            oscillator.pndiode = self

        super(PNdiode,self).__init__(**kwargs)

    def __str__(self):
        return "{}\nDiode:\n Gain = {:e}\n Dark current = {}\nSecondary target:\n {}\nOptics:\n {}\nVoltage-to-Frequency:\n {}"\
                .format(super(PNdiode,self).__str__(),\
                self.Rout,self.darkcurrent,self.secondarytarget,\
                self.optics,self.oscillator)    
    
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

    def _diode_attenuation(self,energy,thickness):
        """Calculate attenuation: 1-exp(-mu.rho.thickness)

        Args:
            energy(num): keV
            thickness: cm

        Returns:
            num or array-like: attenuation
        """
        return (1-np.exp(-self.material.mass_att_coeff(units.magnitude(energy,"keV"))*self.material.density*units.magnitude(thickness,"cm")))
    
    def thicknessfromattenuation(self,energy,attenuation):
        """Calculate thickness: att = 1-exp(-mu.rho.thickness)

        Args:
            energy(num): keV
            attenuation: 

        Returns:
            num or array-like: attenuation
        """
        thickness = -np.log(1-attenuation)/(self.material.mass_att_coeff(units.magnitude(energy,"keV"))*self.material.density)
        self.thickness = ureg.Quantity(thickness,"cm")
    
    def spectral_responsivity(self,energy):
        """Generated current per radiative power

        Args:
            energy(num): keV

        Returns:
            num or array-like: A/W or e/eV
        """
        if self.thickness is None:
            diodeatt = 1
        else:
            diodeatt = self._diode_attenuation(energy,self.thickness)
        
        return diodeatt*(self.ELCHARGE/units.Quantity(self.ehole,"eV")).to("A/W")

    def _chargeperdiodephoton(self,energy):
        """Charge generated per photon absorbed by the diode: spectral responsivity multiplied by the energy

        Args:
            energy(num): keV

        Returns:
            num or array-like: C/ph
        """
        return (self.spectral_responsivity(energy)*units.Quantity(energy,"keV")).to("coulomb")
    
    def _secondarytarget(self,energy,weights):
        """
        Args:
            energy(array): source energies in keV (dims: nE)
            weights(array): source energies in keV (dims: nE)
            
        Returns:
            energy2(array): secondary energies in keV (dims: nE2)
            Ts(array): transmission of source energies (dims: nE)
            wY(array): rates of secondary lines (dims: nE x nE2)
            
        """
        if self.secondarytarget is None:
            Ts = np.ones_like(energy)
            wY = np.ones_like(energy)
        else:
            nE = len(energy)
            
            # transmission of the source spectrum (target and beam filters)
            Ts = self.secondarytarget.transmission(energy)*self.secondarytarget.geometry.filter_transmission(energy,source=True)
            
            # spectrum generated from the target
            spectrum = self.secondarytarget.xrayspectrum(energy,weights=weights) 

            # Extract energies and rates
            geomkwargs = self.secondarytarget.geometry.xrayspectrumkwargs()
            energy,wY = zip(*list(spectrum.spectrum(**geomkwargs)))

            energy = np.asarray(energy)
            nE2 = len(energy)
            wY = np.asarray(wY).reshape((nE,nE2))

        return energy,Ts,wY
        
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
        
    def _chargepersamplephoton(self,energy,weights=None):
        """Charge generated per photon reaching the sample

        Args:
            energy(num|array): source energies in keV (dims: nE)
            weights(num|array): source line weights (dims: nE)

        Returns:
            num: C/ph
        """
        
        # Parse input (remove units and make array)
        energy,func = instance.asarrayf(units.magnitude(energy,"keV"))
        if weights is None:
            weights = np.ones_like(energy)
        else:
            weights,func = instance.asarrayf(weights)
        weights /= weights.sum() # make sure those are fractions
        
        energy2,Ts,wY = self._secondarytarget(energy,weights) # nE2, nE, nE x nE2
        To = self._transmission_optics(energy) # nE
        chargeperdiodephoton = self._chargeperdiodephoton(energy2) # nE2
        
        # I(A)  = Is(ph/s) . Cs(C) + D(A)
        # Cs(C) = SUM_i [SUM_j[wi.Yij.Cj]] / SUM_k [wk.Tks.Tko]
        #
        # wi: source line fraction
        # Cj (C/ph): charge per photon hitting the diode
        # Cs (C/ph): charge per photon hitting the sample
        #

        # Sum over secondary lines: SUM_j[wi.Yij.Cj]
        wY = np.sum(wY*chargeperdiodephoton[np.newaxis,:],axis=1) # nE

        # Sum over source lines: 
        return np.sum(wY)/np.sum(weights*Ts*To)
    
    def samplelineweights(self,energy,weights):
        """Source weights after passing through the diode+optics
        
        Args:
            energy(num|array): source energies in keV (dims: nE)
            weights(num|array): source line weights (dims: nE)

        Returns:
            weights(num|array): source line weights at the sample position
        """
        
        # Parse input (remove units and make array)
        energy,func = instance.asarrayf(units.magnitude(energy,"keV"))
        if weights is None:
            weights = np.ones_like(energy)
        else:
            weights,func = instance.asarrayf(weights)
        weights /= weights.sum() # make sure those are fractions
    
        if self.secondarytarget is None:
            Ts = 1
        else:
            Ts = self.secondarytarget.transmission(energy)
        To = self._transmission_optics(energy)

        # wis = Iis/sum_i[Iis] = wi.Tis.Tio / sum_i[wi.Tis.Tio]
        weights *= Ts*To
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
            spectrocrunch.math.linop: slope (C/ph), intercept(C/s)
        """
        return linop(self._chargepersamplephoton(energy),self._dark())

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
        
    def fluxtocurrent(self,energy,flux):
        """
        Args:
            flux(num or array-like): ph/s
            energy(num): keV

        Returns:
            num or array-like: current (A)
        """

        op = self.op_fluxtocurrent(units.Quantity(energy,"keV"))
        return op(units.Quantity(flux,"Hertz")).to("ampere")

    def currenttoflux(self,energy,current):
        """
        Args:
            flux(num or array-like): A
            energy(num): keV

        Returns:
            num or array-like: flux (ph/s)
        """

        op = self.op_currenttoflux(units.Quantity(energy,"keV"))
        return op(units.Quantity(current,"ampere")).to("hertz")

    def cpstoflux(self,energy,cps):
        """
        Args:
            cps(num or array-like): cts/s
            energy(num): keV

        Returns:
            num or array-like: flux (ph/s)
        """

        op = self.op_cpstoflux(units.Quantity(energy,"keV"))
        return op(units.Quantity(cps,"hertz")).to("hertz")

    def fluxtocps(self,energy,flux):
        """
        Args:
            flux(num or array-like): ph/s
            energy(num): keV

        Returns:
            num or array-like: cps (cts/s)
        """

        op = self.op_fluxtocps(units.Quantity(energy,"keV"))
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

    def xrfnormop(self,energy,time,reference):
        """Operator to convert the raw diode signal to a flux normalizing signal.
           Usage: Inorm = I/op(iodet)
        
            XRF flux normalization:
            
            Measured XRF:
            Ixrf(cts) = F(cps).t(s).cxrf
            
            Desired XRF signal:
            Ixrf(cts)/norm = Fref(cps).t(s).cxrf
            
            Normalization function to be apply on the raw diode signal Idiode
            norm = F/Fref = cpstoflux(Idiode/t)/Fref = op(Idiode)
            op: x-> cpstoflux(x/t)/Fref
            
            In case reference in counts instead of photons/sec
            Fref = round_sig(cpstoflux(Idioderef/t),2)
    
        Args:
            energy(num): keV
            time(num): sec
            reference(num): iodet (counts) or flux (photons/sec) to which the data should be normalized

        Returns:
            op(linop): raw diode conversion operator
            Fref(num): reference flux in photons/s
            time(num): reference time in s
        """

        # Convert from counts to photons/sec
        # op: x-> cpstoflux(x/t)
        op = self.op_cpstoflux(energy)
        op.m /= units.Quantity(time,"s")
        
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
        op.m = units.magnitude(op.m,"dimensionless")
        op.b = units.magnitude(op.b,"dimensionless")
        
        return op,Fref.to("hertz").magnitude,units.magnitude(time,"s")
        
    def gainfromcps(self,energy,cps,fluxest):
        """Estimate the gain, assuming it is 10^x V/A
        
        Args:
            energy(num): keV
            cps(num): 
            fluxest(num): estimated flux 

        Returns:
            None
        """
        fluxcalc = self.cpstoflux(energy,cps)
        r = np.log10(fluxcalc/units.magnitude(units.Quantity(fluxest,"hertz"),"dimensionless"))
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
        
    def setopticstransmission(self,energy,transmission):
        """Set optics transmission for a particular energy
        
        Args:
            energy(num): keV
            transmission(num)

        Returns:
            None
        """
        self.optics.set_transmission(units.magnitude(energy,"keV"),transmission)
        
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
        self.finterpol = interpolate.interp1d(energy,response)

    def spectral_responsivity(self,energy):
        """Return spectral responsivity

        Args:
            energy(num or array-like): keV

        Returns:
            num or array-like: A/W
        """
        if self.model:
            r = self._diode_attenuation(energy,self.thickness)/self.ehole*self.ELCHARGE
        else:
            try:
                r = self.finterpol(energy)
            except:
                r = self._diode_attenuation(energy,self.thickness)/self.ehole*self.ELCHARGE
         
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
        #energy = ureg.Quantity(np.linspace(min(energy)*0.9,max(energy)*1.1,100),'keV')
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
        #energy = ureg.Quantity(np.linspace(min(energy)*0.9,max(energy)*1.1,100),'keV')
        #plt.plot(energy,self._spectral_responsivity_ratio(energy))
        #plt.show()
        

class NonCalibratedPNdiode(PNdiode):

    def __init__(self, **kwargs):
        super(NonCalibratedPNdiode,self).__init__(**kwargs)

    def calibrate(self,cps,sampleflux,energy,caliboption=None):
        """Calibrate with another diode measuring the flux at the sample position

        Args:
            cps(array): count rate measured by this diode (cts/s)
            sampleflux(array): flux measured at the sample position (ph/s)
            energy(num): keV
            caliboption(Optional(string)): "optics" (default), "target" or "diode"
            
        Returns:
            None
        """
        energy = units.Quantity(energy,"keV")
            
        current = self.cpstocurrent(units.Quantity(cps,"hertz"))

        # I(A) = Is(ph/s).slope + intercept
        # I(A) = Is(ph/s).chargepersamplephoton/cor + D(A)
        # chargepersamplephoton = Y/(Ts.To).E(eV/ph).1(e).(1-T)/Ehole(eV)
        
        slope,intercept = linfit(units.magnitude(sampleflux,"hertz"),units.magnitude(current,"ampere"))
        slope = units.Quantity(slope,"ampere/hertz")
        intercept = units.Quantity(intercept,"ampere")
        
        # Adapt dark current:
        self.setdark(intercept)
        
        # Adapt diode thickness, solid angle or transmission
        cor = units.magnitude(self._chargepersamplephoton(energy)/slope,"dimensionless")

        if caliboption=="target":
            # Correct the secondary target yield by changing the diode solid angle
            # Y' = Y/cor
            sa = self.secondarytarget.solidangle/cor
            if sa<0 or sa>(2*np.pi):
                    raise RuntimeError("Diode solid angle of 4*pi*{} srad is not valid (possible wrong parameters: optics transmission, diode thickness)".format(sa/(4*np.pi)))
            self.solidangle = sa
        elif caliboption=="diode":
            # Correct the diode thickness
            # (1-T')/To' = (1-T)/(To.cor)
            # diode transmission (function of thickness): 0 < T' <= 1
            attenuation = self._diode_attenuation(energy,self.thickness)/cor
            if attenuation<0 or attenuation>1:
                raise RuntimeError("Diode attenuation of {} % is not valid (possible wrong parameters: optics transmission, diode solid angle)".format(attenuation*100))
            self.thicknessfromattenuation(energy,attenuation)
        else:
            # Correct the optics transmission
            # (1-T')/To' = (1-T)/(To.cor)
            # optics transmission: 0 < To' <= 1
            if self.optics is None:
                raise RuntimeError("Optics needed for calibration, otherwise specify the \"caliboption\" keyword")
            To = self.optics.transmission(units.magnitude(energy,"keV"))*cor
            if To<0 or To>1:
                raise RuntimeError("Optics transmission of {} % is not valid (possible wrong parameters: diode thickness, diode solid angle)".format(To*100))
            self.setopticstransmission(energy,To)


class SXM_PTB(AbsolutePNdiode):
    aliases = ["ptb"]
    
    def __init__(self,**kwargs):
        attenuators = {}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":np.inf}
        kwargs["attenuators"] = attenuators
        
        ptb = np.loadtxt(resource_filename('id21/ptb.dat'))
        energy = ureg.Quantity(ptb[:,0],"keV")
        response = ureg.Quantity(ptb[:,1],"milliampere/watt")

        super(SXM_PTB,self).__init__(\
                        Rout=ureg.Quantity(1e5,"volt/ampere"),\
                        darkcurrent=ureg.Quantity(0,"ampere"),\
                        energy=energy,response=response,**kwargs)

class SXM_IDET(CalibratedPNdiode):
    # International Radiation Detectors (IRD), AXUV...
    # Keithley K428 (10V max analog output)
    # NOVA N101VTF voltage-to-frequency converter (Fmax=1e6, F0=0Hz)
    # P201 counter board
    
    aliases = ["idet"]
    
    def __init__(self,**kwargs):
        attenuators = {}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":np.inf}
        kwargs["attenuators"] = attenuators
        
        absdiode = SXM_PTB(model=True)

        vtof = Oscillator(Fmax=ureg.Quantity(1e6,"Hz"),\
                        F0=ureg.Quantity(0,"Hz"),\
                        Vmax=ureg.Quantity(10,"volt"))

        ird = np.loadtxt(resource_filename('id21/ird.dat'))
        energy = ureg.Quantity(ird[:-1,0],"keV")
        responseratio = ird[:-1,1]

        super(SXM_IDET,self).__init__(\
                        Rout=ureg.Quantity(1e5,"volt/ampere"),\
                        darkcurrent=ureg.Quantity(0,"ampere"),\
                        energy=energy,reponseratio=responseratio,\
                        absdiode=absdiode,oscillator=vtof,**kwargs)
                        
class SXM_IDET_TEST(AbsolutePNdiode):

    def __init__(self,**kwargs):
        # This is a test: assumes sxmidet responds like an absolute diode (which is not the case)
        attenuators = {}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":np.inf}
        kwargs["attenuators"] = attenuators

        absdiode = sxmptb(model=model)

        vtof = Oscillator(Fmax=ureg.Quantity(1e6,"Hz"),\
                        F0=ureg.Quantity(0,"Hz"),\
                        Vmax=ureg.Quantity(10,"volt"))

        ird = np.loadtxt(resource_filename('id21/ird.dat'))
        energy = ureg.Quantity(ird[:-1,0],"keV")
        responseratio = ird[:-1,1]
        response = responseratio*absdiode.spectral_responsivity(energy)
    
        super(SXM_IDET_TEST,self).__init__(\
                        Rout=ureg.Quantity(1e5,"volt/ampere"),\
                        darkcurrent=ureg.Quantity(0,"ampere"),\
                        energy=energy,response=response,\
                        oscillator=vtof,**kwargs)

class SXM_IODET1(NonCalibratedPNdiode):
    # International Radiation Detectors (IRD), AXUV-PS1-S
    # Keithley K428 (10V max analog output)
    # NOVA N101VTF voltage-to-frequency converter (Fmax=1e6, F0=0Hz)
    # P201 counter board

    aliases = ["iodet1"]

    def __init__(self,optics=True):
        kwargs = {}
        
        attenuators = {}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":14e-4}
        kwargs["attenuators"] = attenuators
        
        source = xraysources.factory("synchrotron")
        geometry = diodegeometries.factory("Geometry",anglein=90,angleout=45,source=source,detector=self)
        kwargs["solidangle"] = 0.4*4*np.pi
        kwargs["ehole"] = constants.eholepair_si()
        
        window = compoundfromname.compoundfromname("silicon nitride")
        coating = element.Element('Si')
        secondarytarget = multilayer.Multilayer(material=[coating,window],thickness=[0.5e-4,0.5e-4],geometry=geometry)

        if optics:
            optics = KB.KB()
        else:
            optics = None
            
        vtof = Oscillator(Fmax=ureg.Quantity(1e6,"Hz"),\
                            F0=ureg.Quantity(0,"Hz"),\
                            Vmax=ureg.Quantity(10,"volt"))
                        
        super(SXM_IODET1,self).__init__(\
                            Rout=ureg.Quantity(1e5,"volt/ampere"),\
                            darkcurrent=ureg.Quantity(0,"ampere"),\
                            oscillator=vtof,\
                            secondarytarget=secondarytarget,\
                            optics=optics,\
                            **kwargs)
                        
class SXM_IODET2(NonCalibratedPNdiode):
    # International Radiation Detectors (IRD), AXUV-PS1-S
    # Keithley K428 (10V max analog output)
    # NOVA N101VTF voltage-to-frequency converter (Fmax=1e6, Vmax=0Hz)
    # P201 counter board
    
    aliases = ["iodet2"]

    def __init__(self,optics=True):
        kwargs = {}
        
        attenuators = {}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":14e-4}
        kwargs["attenuators"] = attenuators
        
        source = xraysources.factory("synchrotron")
        geometry = diodegeometries.factory("Geometry",anglein=90,angleout=45,source=source,detector=self)
        kwargs["activearea"] = 0.4*4*np.pi
        kwargs["ehole"] = constants.eholepair_si()
        
        window = compoundfromname.compoundfromname("silicon nitride")
        secondarytarget = multilayer.Multilayer(material=[window],thickness=[0.5e-4],geometry=geometry)
        
        if optics:
            optics = KB.KB()
        else:
            optics = None
            
        vtof = Oscillator(Fmax=ureg.Quantity(1e6,"Hz"),\
                            F0=ureg.Quantity(0,"Hz"),\
                            Vmax=ureg.Quantity(10,"volt"))
                        
        super(SXM_IODET2,self).__init__(\
                            Rout=ureg.Quantity(1e5,"volt/ampere"),\
                            darkcurrent=ureg.Quantity(0,"ampere"),\
                            oscillator=vtof,\
                            secondarytarget=secondarytarget,\
                            optics=optics,\
                            **kwargs)

class XRD_IDET(NonCalibratedPNdiode):
    aliases = ["pico1"]
    
    def __init__(self,**kwargs):
        attenuators = {}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":0.1}
        kwargs["attenuators"] = attenuators
        
        kwargs["ehole"] = constants.eholepair_si()
        
        
        super(XRD_IDET,self).__init__(\
                                Rout=ureg.Quantity(10,"volt")/ureg.Quantity(2.1e-6,"ampere"),\
                                darkcurrent=ureg.Quantity(0,"ampere"),
                                **kwargs)

factory = PNdiode.factory


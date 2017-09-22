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

from __future__ import division

from . import areadetectors
from . import scintillators
from . import lenses
from . import noisepropagation
from ..common.instance import isarray

import numpy as np

import uncertainties

import itertools

import matplotlib.pyplot as plt

def transmission(N,N0,D=0,D0=0,\
                tframe_data=1,nframe_data=1,\
                tframe_flat=1,nframe_flat=1,\
                nframe_dark=1):
    num = N - D/nframe_dark*nframe_data
    num /= nframe_data*tframe_data
    denom = N0 - D0/nframe_dark*nframe_flat
    denom /= nframe_flat*tframe_flat
    return num/denom
    
def absorbance(_transmission):
    ret = np.empty_like(_transmission)
    ind = _transmission>0
    ret[ind] = -uncertainties.unumpy.log(_transmission[ind])
    ret[~ind] = noisepropagation.randomvariable(-np.inf,np.inf)
    return ret

class id21_ffsetup(object):

    def __init__(self,composition=None,scint="LSO",lens="x10"):
        self.composition = composition

        if scint=="LSO":
            self.oscint = scintillators.factory("LSO ID21",thickness=10)
        else:
            self.oscint = scintillators.factory("GGG ID21",thickness=13)

        if lens=="x10":
            self.olens = lenses.factory("mitutoyoid21_10x")
        else:
            self.olens = lenses.factory("mitutoyoid21_10x")

        self.odet = areadetectors.factory("pcoedge55")
    
    def propagate(self,E,N,tframe,nframe,samplein=False):
        if isarray(E):
            E = np.asarray(E)
        
        if samplein and self.composition is not None:
            N = self.composition.propagate(N,E)
        N = self.oscint.propagate(N,E)
        E = self.oscint.visspectrum
        
        if isarray(N):
            s = N.shape
            N = N.flatten()
        N = self.olens.propagate(N,E,nrefrac=self.oscint.get_nrefrac())
        N = self.odet.propagate(N,E,tframe=tframe,nframe=nframe)
        
        if isarray(N):
            N = N.reshape(s)
        return N
    
    def measurement(self,flux,energy,\
                    tframe_data=None,nframe_data=None,\
                    tframe_flat=None,nframe_flat=None,\
                    nframe_dark=None):
        """ID21 fullfield noise propagation

        Args:
            flux (num|array): incomming flux (ph/sec)
            energy (num|array): associated energy (keV)
            sample (spectrocrunch.simulation.materials.Material): sample
            tframe_data(num): time per frame (sec)
            nframe_data(num): number of data frames
            tframe_flat(num): time per frame (sec)
            nframe_flat(num): number of flat frames
            nframe_dark(num): number of dark frames
            scint(Optional(str)):
            lens(Optional(str)):

        Returns:
            4-tuple(uncertainties.unumpy.uarray): detector signal in DU (with and w.o. sample, with and w.o. exposure)
        """
    
        # With sample
        N = flux*tframe_data
        N = noisepropagation.poisson(N)
        N = self.propagate(energy,N,tframe_data,nframe_data,samplein=True)
        
        # Without sample
        N0 = flux*tframe_flat
        N0 = noisepropagation.poisson(N0)
        N0 = self.propagate(energy,N0,tframe_flat,nframe_flat,samplein=False)
        
        # Without beam
        if nframe_dark==0:
            D = noisepropagation.randomvariable(0,0)
            D0 = noisepropagation.randomvariable(0,0)
        else:
            D = self.odet.propagate(noisepropagation.randomvariable(0,0),energy,tframe=tframe_data,nframe=nframe_dark)
            D0 = self.odet.propagate(noisepropagation.randomvariable(0,0),energy,tframe=tframe_flat,nframe=nframe_dark)
    
        return N,N0,D,D0
    
    def fluxquant(self,energy,tframe_flat=None,nframe_flat=None,nframe_dark=None):
        """ From a dark subtract flat field image (DU/sec) we can calculate the flux as follows:
                flux = img/fluxquant
        
        Args:
            energy(num):
        
        Returns:
            num: fluxquant (DU/ph)
        """
        
        # Without sample
        N0 = noisepropagation.randomvariable(1,0)
        N0 = self.propagate(energy,N0,tframe_flat,nframe_flat,samplein=False)
        D0 = self.odet.propagate(noisepropagation.randomvariable(0,0),energy,tframe=tframe_flat,nframe=nframe_dark)
        ret = N0/nframe_flat
        if nframe_dark!=0:
            ret -= D0/nframe_dark
        return noisepropagation.E(ret)
    
    @staticmethod
    def getnframes(totaltime,frametime,fracflat):
        n = int(round(totaltime/frametime))
        nflat = max(int(round(fracflat*n/2.)),1)
        nflat *= 2 # before and after
        ndata = max(n - nflat,1)

        return ndata,nflat

    @staticmethod
    def getrealtime(totaltime,frametime,fracflat):
        ndata,nflat = self.getnframes(totaltime,frametime,fracflat)
        n = ndata+nflat
        overhead = 6.50305 + 0.0131498*n
        return frametime*n + overhead

    def xanes(self,flux,energy,**kwarg):

        N,N0,D,D0 = self.measurement(flux,energy,**kwarg)
                                
        T = transmission(N,N0,D=D,D0=D0,**kwarg)

        XAS = absorbance(T)

        signal = noisepropagation.E(XAS)
        noise = noisepropagation.S(XAS)
        return signal,noise
    
    def costfunc(self,flux,energy,**kwargs):
        signal,noise = self.xanes(flux,energy,**kwargs) 
        return np.max(noise)/(signal[-1]-signal[0])
    
    def __str__(self):
        return str(self.composition)
    
    def plotxanesnoise(self,flux,energy,**kwargs):
        signal,noise = self.xanes(flux,energy,**kwargs)
        signal = np.random.normal(signal,noise)
        plt.plot(energy,signal)
        plt.xlabel("Energy (keV)")
        plt.ylabel("N/S (%)")
       
    def plotxanesNSR(self,flux,energy,**kwargs):
        signal,noise = self.xanes(flux,energy,**kwargs)
        plt.plot(energy,noise/signal*100)
        plt.xlabel("Energy (keV)")
        plt.ylabel("N/S (%)")
   
    def plotxanes(self,flux,energy,**kwargs):
        signal,_ = self.xanes(flux,energy,**kwargs)
        plt.plot(energy,signal)
        plt.xlabel("Energy (keV)")
        plt.ylabel("Absorbance")

def id21_transmissionnoise(flux,energy,time,iodetgain,idetgain,iodet="1"):
    """ID21 transmission noise propagation
    """
    pass

def id21_xrfnoise(flux,energy,time):
    """ID21 XRF noise propagation
    """
    pass



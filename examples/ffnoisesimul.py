# -*- coding: utf-8 -*-

import os, sys
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spectrocrunch.materials.compoundfromformula import compoundfromformula
from spectrocrunch.materials.compoundfromname import compoundfromname
from spectrocrunch.materials.mixture import mixture
from spectrocrunch.materials.types import fractionType

from spectrocrunch.simulation import calcnoise
from spectrocrunch.simulation import materials
from spectrocrunch.simulation import noisepropagation

import numpy as np
import scipy.optimize

import matplotlib.pyplot as plt

class sample(object):

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

    def xanes(self,I0,energy,totaltime=None,frametime=None,fracflat=None,ndark=None):
        ndata,nflat = self.getnframes(totaltime,frametime,fracflat)

        energy = np.asarray(energy)

        N,N0,D,D0 = calcnoise.id21_ffnoise(I0,energy,self.composition,\
                                tframe_data=frametime,nframe_data=ndata,\
                                tframe_flat=frametime,nframe_flat=nflat,\
                                nframe_dark=ndark)
                                
        T = calcnoise.transmission(N,N0,D=D,D0=D0,\
                tframe_data=frametime,nframe_data=ndata,\
                tframe_flat=frametime,nframe_flat=nflat,\
                nframe_dark=ndark)

        XAS = calcnoise.absorbance(T)

        signal = noisepropagation.E(XAS)
        noise = noisepropagation.S(XAS)
        return signal,noise
    
    def costfunc(self,I0,energy,**kwargs):
        signal,noise = self.xanes(I0,energy,**kwargs) 
        
        #return np.max(noise/signal*100)

        return np.mean(noise)/(signal[-1]-signal[0])
        
    
    def __str__(self):
        return str(self.composition)
        
    def plotxanesnoise(self,I0,energy,**kwargs):
        signal,noise = self.xanes(I0,energy,**kwargs)
        plt.plot(energy,noise/signal*100)
        plt.xlabel("Energy (keV)")
        plt.ylabel("N/S (%)")
   
    def plotxanes(self,I0,energy,**kwargs):
        signal,_ = self.xanes(I0,energy,**kwargs)
        plt.plot(energy,signal)
        plt.xlabel("Energy (keV)")
        plt.ylabel("Absorbance")
             
class sample_hg115(sample):

    def __init__(self,wpigment=10,paintthickness=10):
        binder = compoundfromname("linseed oil")
        pigment = compoundfromname("verdigris")
        paint = mixture([binder,pigment],[1-wpigment/100.,wpigment/100.],fractionType.weight)
        ultralene = compoundfromname("ultralene")
        sfreetape = compoundfromname("sulfur-free tape")

        #ultralene = compoundfromname("vacuum")
        #sfreetape = compoundfromname("vacuum")
        
        m = [ultralene,paint,sfreetape]
        thickness = [4,paintthickness,10]
        
        #m = [compoundfromname("vacuum"),compoundfromname("vacuum"),compoundfromname("vacuum")]

        self.composition = materials.factory("Multilayer",material=m,thickness=thickness,anglein=0,angleout=0,azimuth=0)
        self.paintindex = 1

    def set_wpigment(self,wpigment):
        w = self.composition.material[self.paintindex].weightfractions()
        w["verdigris"] = wpigment/100.
        w["linseed oil"] = 1-wpigment/100.
        self.composition.material[self.paintindex].change_fractions(w,fractionType.weight)
        
    def get_wpigment(self):
        return self.composition.material[self.paintindex].weightfractions()["verdigris"]*100
    
    def set_paintthickness(self,paintthickness):   
        self.composition.thickness[self.paintindex] = paintthickness

    def get_paintthickness(self):
        return self.composition.thickness[self.paintindex]
        
    def optimize_thickness(self,I0,energy,**kwargs):
        def costfunc(paintthickness):
            self.set_paintthickness(paintthickness[0])
            c = self.costfunc(I0,energy,**kwargs)
            return c
        
        guess = self.get_paintthickness()
        
        result = scipy.optimize.least_squares(costfunc, guess, gtol=1e-015, ftol=1e-015)
        print result.message
        
        return result.x[0],result.success

    def optimize_wpigment(self,I0,energy,**kwargs):
        def costfunc(wpigment):
            self.set_wpigment(wpigment[0])
            c = self.costfunc(I0,energy,**kwargs)
            return c
        
        guess = self.get_wpigment()
        
        result = scipy.optimize.least_squares(costfunc, guess, bounds=([0,100]), gtol=1e-015, ftol=1e-015)
        print result.message
        
        return result.x[0],result.success
        
    def optimize_thickness_plot(self,I0,energy,**kwargs):
        thickness = self.get_paintthickness()
        
        t = np.linspace(max(thickness-100,0),thickness+100,50)
        r = np.zeros(len(t))
        for i,paintthickness in enumerate(t):
            self.set_paintthickness(paintthickness)
            r[i] = self.costfunc(I0,energy,**kwargs)

        self.set_paintthickness(thickness)

        plt.plot(t,1/r,'-o',label="{} %".format(self.get_wpigment()))
        plt.xlabel("thickness ($\mu$m)")
        plt.ylabel("Jump-to-noise")
        
    def optimize_wpigment_plot(self,I0,energy,**kwargs):
        w = self.get_wpigment()
    
        t = np.linspace(0,20,50)
        r = np.zeros(len(t))
        for i,wpigment in enumerate(t):
            self.set_wpigment(wpigment)
            r[i] = self.costfunc(I0,energy,**kwargs)

        self.set_wpigment(w)
        
        plt.plot(t,1/r,'-o',label="{} $\mu$m".format(self.get_paintthickness()))
        plt.xlabel("Verdigris (%)")
        plt.ylabel("Jump-to-noise")
            
    def optimize(self,I0,energy,**kwargs):

        def costfunc(p):
            self.set_wpigment(p[0])
            self.set_paintthickness(p[1])
            return self.costfunc(I0,energy,**kwargs)
        
        guess = (self.get_wpigment(),self.get_paintthickness())
        
        result = scipy.optimize.least_squares(costfunc, guess, bounds=([0,0],[100,1e6]), gtol=1e-015)
        print result.message
        
        return result.x,result.success 
        
def hg115_ff():
    sample = sample_hg115()

    I0 = 1e6
    energy = np.linspace(8.9,9.3,100)
    totaltime = 70
    frametime = 0.07
    fracflat = 1/3.
    ndark = 30
    
    kwargs = {"totaltime":totaltime,\
                "frametime":frametime,\
                "fracflat":fracflat,
                "ndark":ndark}

    opt = 1

    energyopt = [8.97,9]
    if opt==0:
        sample.set_wpigment(10)
        t,s = sample.optimize_thickness(I0,energyopt,**kwargs)
        sample.set_paintthickness(t)
    elif opt==1:
        sample.set_paintthickness(20)
        w,s = sample.optimize_wpigment(I0,energyopt,**kwargs)
        sample.set_wpigment(w)
    else:
        wt,s = sample.optimize(I0,energy,**kwargs)
        sample.set_wpigment(wt[0])
        sample.set_paintthickness(wt[1])
    
    print "Thickness = {} μm".format(sample.get_paintthickness())
    print "Verdigris = {} wt%".format(sample.get_wpigment())
    print "Jump to noise = {}".format(1/sample.costfunc(I0,energyopt,**kwargs))
    print ""
        
        
    plt.figure()
    for thickness in [10,15,20]:
        sample.set_paintthickness(thickness)
        sample.optimize_wpigment_plot(I0,energy,**kwargs)
    plt.legend(loc='best')
    plt.show()
    
    exit()
    
    sample.optimize_thickness_plot(I0,energy,**kwargs)
    
    
    sample.optimize_wpigment_plot(I0,energy,**kwargs)
                                      
    plt.figure()
    sample.plotxanes(I0,energy,**kwargs)
      
    plt.figure()
    sample.plotxanesnoise(I0,energy,**kwargs)    
          
    plt.show()

def hg115_xrd():
    sample = sample_hg115()
    energy = 8.5
    sample.set_wpigment(100)
    
    r = np.linspace(10,20,50)
    n = [None]*len(r)
    for i,t in enumerate(r):
        sample.set_paintthickness(t)
        n[i] = noisepropagation.E(sample.composition.propagate(noisepropagation.poisson(1e7),energy,interaction=materials.interactionType.elastic))

    print n[-1]/n[0]
    plt.plot(r,n)
    plt.show()

if __name__ == '__main__':
    hg115_ff()


#        I0 = 1e5
#        energy = np.linspace(3,5,100)
#        tframe = 0.07
#        nframe = 100
#        ndark = 30
#tframe_data=tframe,nframe_data=nframe,\
#                    tframe_flat=tframe,nframe_flat=nframe,\
#                    nframe_dark=ndark

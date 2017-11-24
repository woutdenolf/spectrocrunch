# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
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

import spectrocrunch.io.spec as spec
import spectrocrunch.detectors.diode as diode
import matplotlib.pyplot as plt

class FluxMonitor(object):
    
    def __init__(self):
        self.idet = diode.factory("idet",model=True)
        self.iodet = None
        
    def parsespec(self,specfile,specnr):
        # Parse spec file
        fspec = spec.spec(specfile)
        
        self.ioz,self.istopz,self.energy,self.zpz = fspec.getmotorvalues(specnr,["diodeIoZ","istopz","Energy MONO",'zpz'])

        data = fspec.getdata2(specnr,["Seconds","Photons","idet","iodet"])
        seconds = data[:,0]
        self.fluxest = data[:,1]/seconds
        self.idetcps = data[:,2]/seconds
        self.iodetcps = data[:,3]/seconds
        
    def setiodet(self):
        if abs(self.ioz-7)<abs(self.ioz-23):
            name = "iodet1"
        elif abs(self.istopz+20)<abs(self.istopz+1.3):
            name = "iodet2"
        else:
            raise RuntimeError("No diode in the beam")

        if self.iodet is None:
            optics = (abs(self.zpz-7)<abs(self.zpz-6.5))
            
            self.iodet = diode.factory(name,optics=optics)
            self.iodet.name = name
            #self.iodet.secondarytarget.geometry.solidangle = 4*np.pi*0.2
        else:
            if self.iodet.name!=name:
                raise RuntimeError("Wrong diode in the beam")

    def darkiodet(self):
        self.iodet.darkfromcps(self.iodetcps)

    def darkidet(self):
        self.idet.darkfromcps(self.idetcps)

    def checkdark(self):
        f, (ax1, ax2) = plt.subplots(1, 2)
        self._checkdark(self.idet,self.idetcps,"idet",ax1)
        self._checkdark(self.iodet,self.iodetcps,"iodet",ax2)
        plt.tight_layout()
        
    def _checkdark(self,d,cps,name,ax):
        ax.plot(cps,'o',label='spec')
        ax.axhline(y=d.fluxtocps(self.energy,0).magnitude,label='calc')
        ax.set_xlabel("points")
        ax.set_ylabel("{} (cts/s)".format(name))
        ax.set_title("{:.0}".format(d.Rout))
        ax.legend(loc="best")
        
    def calibiodet(self):
        if True:
            # Reset the idet gain based on the flux given by spec (fluxest)
            #self.idet.gainfromcps(self.energy,self.idetcps,self.fluxest)
            
            flux = self.idet.cpstoflux(self.energy,self.idetcps)
        else:
            flux = self.fluxest
        
        self.iodet.calibrate(self.iodetcps,flux,self.energy)
        
    def checkflux(self):
        f, (ax1, ax2) = plt.subplots(1, 2)
        
        flux = self.idet.cpstoflux(self.energy,self.idetcps)
        
        ax1.plot(self.idetcps,self.fluxest,'o',label='spec')
        ax1.plot(self.idetcps,flux,label='calc idet')
        ax1.set_xlabel("idet (cts/s)")
        ax1.set_ylabel("flux (ph/s)")
        ax1.set_title("{:.0e}".format(self.idet.Rout))
        ax1.legend(loc="best")
        
        ax2.plot(self.iodetcps,flux,'o',label='calc idet')
        ax2.plot(self.iodetcps,self.iodet.cpstoflux(self.energy,self.iodetcps).magnitude,label='calc iodet')
        ax2.set_xlabel("iodet (cts/s)")
        ax2.set_ylabel("flux (ph/s)")
        ax2.set_title("{:.0e}".format(self.iodet.Rout))
        ax2.legend(loc="best")
    
        plt.tight_layout()
    
    def initdiodes(self,gainiodet,gainidet):
        self.iodet.setgain(gainiodet)
        self.idet.setgain(gainidet)
        
    def initspec(self,specfile,specnr,gainiodet,gainidet):
        self.parsespec(specfile,specnr)
        self.setiodet()
        self.initdiodes(gainiodet,gainidet)

    def dark(self,specfile,specnr,gainiodet,gainidet,plot=False):
        self.initspec(specfile,specnr,gainiodet,gainidet)
        
        self.darkiodet()
        self.darkidet()
        
        if plot:
            self.checkdark()
            plt.show()
 
    def calib(self,specfile,specnr,gainiodet,gainidet,plot=False):
        self.initspec(specfile,specnr,gainiodet,gainidet)

        self.calibiodet()
        
        if plot:
            self.checkflux()
            plt.show()

    def manualcalib(self,energy,transmission,gainiodet,gainidet):
        self.initdiodes(gainiodet,gainidet)
        self.iodet.optics.set_transmission(energy,transmission)
    
    def xrfnormop(self,energy,time,reference):
        return self.iodet.xrfnormop(energy,time,reference)
        

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
import logging

import spectrocrunch.io.spec as spec
import spectrocrunch.detectors.diode as diode
import matplotlib.pyplot as plt

from ..common import units

logger = logging.getLogger(__name__)

class FluxMonitor(object):
    
    def __init__(self,iodetname=None,focussed=None):
        self.idet = diode.factory("idet",model=True)
        self.iodet = None
        self.iodetname = iodetname
        self.focussed = focussed
        self.reference = units.Quantity(1e10,"hertz")
        self.time = units.Quantity(0.1,"s")
        self.setiodet()
        
    def __getattr__(self,attr):
        return getattr(self.iodet,attr)
    
    def __str__(self):
        return "Default flux reference:\n {:~}\nDefault exposure time:\n {:~}\n{}".format(self.reference,self.time,self.iodet.geometry)

    def parsespec(self,specfile,specnr):
        fspec = spec.spec(specfile)
        ioz,istopz,energy,zpz = fspec.getmotorvalues(specnr,["diodeIoZ","istopz","Energy MONO",'zpz'])

        self.checkiodet(ioz,istopz,zpz)

        data = fspec.getdata2(specnr,["Seconds","Photons","idet","iodet"])
        seconds = data[:,0]
        
        self.energy = energy
        self.fluxest = data[:,1]/seconds
        self.idetcps = data[:,2]/seconds
        self.iodetcps = data[:,3]/seconds

    def checkiodet(self,ioz,istopz,zpz):
        if abs(ioz-7)<abs(ioz-23):
            name = "iodet1"
        elif abs(istopz+20)<abs(istopz+1.3):
            name = "iodet2"
        else:
            raise RuntimeError("No diode in the beam")
        if self.iodetname is None:
            self.iodetname = name
        else:
            if name != self.iodetname:
                msg = "Current I0 diode is \"{}\", not {}".format(self.iodetname,name)
                raise RuntimeError(msg)
        
        focussed = abs(zpz-7)<abs(zpz-6.5)
        if self.focussed is None:
            self.focussed = focussed
        else:
            if focussed != self.focussed:
                if self.focussed:
                    msg = "Optics should be in the beam"
                else:
                    msg = "Optics should not be in the beam"
                raise RuntimeError(msg)
                
        self.setiodet()
        
    def setiodet(self):
        if self.iodet is None and not (self.iodetname is None or self.focussed is None):
            self.iodet = diode.factory(self.iodetname,optics=self.focussed)

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
        
    def calibiodet(self,caliboption=None):
        flux = self.idet.cpstoflux(self.energy,self.idetcps)
        #flux = self.fluxest # this is the flux calculated in spec
        
        self.iodet.calibrate(self.iodetcps,flux,self.energy,caliboption=caliboption)
        
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
    
    def prepareddata(self,specfile=None,specnr=None,gainiodet=None,gainidet=None):
        setiodet = False
        setidet = False
        
        loaddata = specfile is not None and specnr is not None
        if loaddata:
            self.parsespec(specfile,specnr)

        if gainiodet is not None:
            self.iodet.setgain(gainiodet)
            setiodet = True
            
        if gainidet is None:
            # Reset the idet gain based on the flux given by spec (fluxest)
            if loaddata:
                self.idet.gainfromcps(self.energy,self.idetcps,self.fluxest)
                setidet = True
        else:
            self.idet.setgain(gainidet)
            setidet = True
            
        return loaddata,setiodet,setidet

    def dark(self,plot=False,**kwargs):
        """Get the dark current from the data
        """
        loaddata,setiodet,setidet = self.prepareddata(**kwargs)
        self.darkiodet()
        self.darkidet()
        
        if plot:
            self.checkdark()
            plt.show()
 
    def setdark(self,iodetcps,idetcps,**kwargs):
        """Manually specify the dark current
        """
        loaddata,setiodet,setidet = self.prepareddata(**kwargs)
        if iodetcps is not None:
            if not setiodet:
                raise RuntimeError("Gain of iodet should be specified")
            self.iodet.darkfromcps(iodetcps)
        if idetcps is not None:
            if not setidet:
                raise RuntimeError("Gain of idet should be specified")
            self.idet.darkfromcps(idetcps)

    def calib(self,caliboption=None,plot=False,**kwargs):
        """Diode calibration from the data
        """
        loaddata,setiodet,setidet = self.prepareddata(**kwargs)
        if not loaddata:
            raise RuntimeError("Attenuator scan should be specified")
        if not setiodet:
            raise RuntimeError("Gain of iodet should be specified")
        if not setidet:
            raise RuntimeError("Gain of idet should be specified")
        self.calibiodet(caliboption=caliboption)
        
        if plot:
            self.checkflux()
            plt.show()

    def setcalib(self,energy,transmission,**kwargs):
        """Manual diode calibration
        """
        loaddata,setiodet,setidet = self.prepareddata(**kwargs)
        self.iodet.optics.set_transmission(energy,transmission)
    
    def setreferenceflux(self,flux):
        self.reference = units.Quantity(flux,"hertz")
    
    def setreferencecounts(self,cts):
        self.reference = units.Quantity(cts,"dimensionless")
    
    def settime(self,time):
        self.time = units.Quantity(time,"s")

    def xrfnormop(self,energy,time=None,ref=None):
        if time is None:
            time = self.time
        if ref is None:
            ref = self.reference
            
        return self.iodet.xrfnormop(energy,time,ref)
    
    def fluxtocps(self,energy,flux):
        return self.iodet.fluxtocps(energy,flux).magnitude
  
    def cpstoflux(self,energy,cps):
        return self.iodet.cpstoflux(energy,cps).magnitude
        

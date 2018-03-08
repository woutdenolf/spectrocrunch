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

import os
import logging
import numpy as np
import matplotlib.pyplot as plt

from ..io import spec
from ..detectors import diode
from ..common import units
from ..detectors import xrf as xrfdetectors
from ..geometries import xrf as xrfgeometries
from ..sources import xray as xraysources

logger = logging.getLogger(__name__)

class FluxMonitor(object):
    
    def __init__(self,iodetname=None,focussed=None,xrfdetector=None,xrfgeometry=None,simplecalibration=True):
        self.source = xraysources.factory("synchrotron")
    
        self.idet = diode.factory("idet",model=False)
        self.iodet = None
        self.simplecalibration = simplecalibration
        self.iodetname = iodetname
        self.focussed = focussed
        self.reference = units.Quantity(1e10,"hertz")
        self.defaulttime = units.Quantity(0.1,"s")
        self.referencetime = None
        self.setiodet()
        
        detector = xrfdetectors.factory(xrfdetector)
        self.xrfgeometry = xrfgeometries.factory(xrfgeometry,detector=detector,source=self.source)
    
    def setxrfposition(self,value):
        self.xrfgeometry.detectorposition = value
    
    def getxrfdistance(self):
        return self.xrfgeometry.distance
    
    def __getattr__(self,attr):
        return getattr(self.iodet,attr)
    
    def __str__(self):
        referencetime = self.referencetime
        if referencetime is not None:
            referencetime = "{~}".format(referencetime)
        return "IODET:\n{}\n\nXRF:\n{}\n\nDefault flux reference:\n {:~e}\Reference exposure time:\n {:~}\nDefault exposure time:\n {:~}\n".\
        format(self.iodet.geometry,self.xrfgeometry,self.reference,referencetime,self.defaulttime)

    def parsespec(self,specfile,specnr):
        fspec = spec.spec(specfile)
        ioz,istopz,energy,zpz = fspec.getmotorvalues(specnr,["diodeIoZ","istopz","Energy MONO",'zpz'])

        self._checkiodet(ioz,istopz,zpz)

        menergies = fspec.haslabel(specnr,"EnergyM")
        if menergies:
            data = fspec.getdata2(specnr,["Seconds","Photons","idet","iodet","EnergyM"])
        else:
            data = fspec.getdata2(specnr,["Seconds","Photons","idet","iodet"])
        seconds = data[:,0]
        self.fluxest = data[:,1]/seconds
        self.idetcps = data[:,2]/seconds
        self.iodetcps = data[:,3]/seconds

        if menergies:
            self.energy = data[:,4]
        else:
            self.energy = energy

    def _checkiodet(self,ioz,istopz,zpz):
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
                msg = "Expected diode is \"{}\", not \"{}\"".format(self.iodetname,name)
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
            self.iodet = diode.factory(self.iodetname,optics=self.focussed,
                                    source=self.source,simplecalibration=self.simplecalibration)

    def darkiodet(self):
        self.iodet.darkfromcps(self.iodetcps)

    def darkidet(self):
        self.idet.darkfromcps(self.idetcps)

    def checkdark(self,plot=False,out=False):
        if plot:
            f, (ax1, ax2) = plt.subplots(1, 2)
            self._checkdark(self.idet,self.idetcps,"idet",ax1)
            color = next(ax1._get_lines.prop_cycler)['color']
            self._checkdark(self.iodet,self.iodetcps,"iodet",ax2,color=color)
            plt.tight_layout()
        if out:
            print("idet :\n Gain = {:~e}\n Dark = {:~e}".format(self.idet.Rout,self.idet.darkcurrent))
            print("iodet:\n Gain = {:~e}\n Dark = {:~e}".format(self.iodet.Rout,self.iodet.darkcurrent))

    def _checkdark(self,d,cps,name,ax,color=None):
        ax.plot(cps,'o',label='spec',color=color)
        ax.axhline(y=d.fluxtocps(5,0).magnitude,label='calc',color=color)
        ax.set_xlabel("points")
        ax.set_ylabel("{} (cts/s)".format(name))
        ax.set_title("{:.0e}".format(d.Rout))
        ax.legend(loc="best")

    def checkflux(self,plot=False,out=False,fluxmin=0,fluxmax=np.inf,fitinfo={}):
        if plot:
            f, (ax1, ax2) = plt.subplots(1, 2)
            
            flux = self.idet.cpstoflux(self.energy,self.idetcps)
            
            ax1.plot(self.fluxest,self.idetcps,'o',label='spec')
            lines = ax1.plot(flux,self.idetcps,label='idet')
            color1 = lines[0].get_color()
            color2 = next(ax1._get_lines.prop_cycler)['color']
            ax1.set_xlabel("flux (ph/s)")
            ax1.set_ylabel("idet (cts/s)")
            ax1.set_title("{:.0e}".format(self.idet.Rout))
            ax1.legend(loc="best")
            
            xoff = 0.3
            yoff = 0.1
            info = self.idet.fluxcpsinfo(self.energy)
            for k,v in info.items():
                ax1.text(xoff,yoff,"{} = {}".format(k,v),transform = ax1.transAxes)
                yoff -= 0.03
            
            
            x = units.magnitude(flux,"hertz")
            indfit = (x<units.magnitude(fluxmin,"hertz")) | (x>units.magnitude(fluxmax,"hertz"))
            if any(indfit):
                ax2.plot(x[indfit],self.iodetcps[indfit],'o',color="#cccccc")
            indfit = np.logical_not(indfit)
            if any(indfit):
                ax2.plot(x[indfit],self.iodetcps[indfit],'o',label='idet',color=color1)
                
            ax2.plot(self.iodet.cpstoflux(self.energy,self.iodetcps).magnitude,self.iodetcps,label='iodet',color=color2)
            ax2.set_xlabel("flux (ph/s)")
            ax2.set_ylabel("iodet (cts/s)")
            ax2.set_title("{:.0e}".format(self.iodet.Rout))
            ax2.legend(loc="best")
            
            if fitinfo:
                xoff = 0.3
                yoff = 0.1
                for k,v in fitinfo.items():
                    ax2.text(xoff,yoff,"{} = {}".format(k,v),transform = ax2.transAxes)
                    yoff -= 0.03
                
            plt.tight_layout()
        
        if out:
            if self.caliboption=="thickness":
                s = "Thickness = {} um".format(self.thickness*1e4)
            elif self.caliboption=="solidangle":
                s = "Solid angle = 4*pi*{} srad".format(self.iodet.solidangle/(4*np.pi))
            elif self.caliboption=="optics":
                s = "Optics transmission = {} %".format(self.iodet.optics.transmission(self.energy)*100)
                
            print("iodet:\n Gain = {:~e}\n Energy = {} keV\n Dark = {:~e}\n {}".format(self.iodet.Rout,self.energy,self.iodet.darkcurrent,s))
        
    def prepareddata(self,specfile=None,specnr=None,gainiodet=None,gainidet=None):
        setiodet = False
        setidet = False
        
        loaddata = specfile is not None and specnr is not None
        if loaddata:
            self.parsespec(specfile,specnr)

        if gainiodet is None:
            # Reset the iodet gain based on the flux given by spec (fluxest)
            if loaddata:
                self.iodet.gainfromcps(self.energy,self.iodetcps,self.fluxest)
                setiodet = True
        else:
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

    def dark(self,plot=False,nofit=False,out=False,**kwargs):
        """Get the dark current from the data
        """
        loaddata,setiodet,setidet = self.prepareddata(**kwargs)
        if not nofit:
            self.darkiodet()
            self.darkidet()
        
        self.checkdark(plot=plot,out=out)
        if plot:
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

    def calib(self,caliboption="optics",fixdark=False,fluxmin=0,fluxmax=np.inf,plot=False,out=False,nofit=False,**kwargs):
        """Diode calibration from the data
        """
        loaddata,setiodet,setidet = self.prepareddata(**kwargs)
        if not loaddata:
            raise RuntimeError("Attenuator scan should be specified")
        if not setiodet:
            raise RuntimeError("Gain of iodet should be specified")
        if not setidet:
            raise RuntimeError("Gain of idet should be specified")
        if not nofit:
            fitinfo = self.calibiodet(caliboption=caliboption,fixdark=fixdark,fluxmin=fluxmin,fluxmax=fluxmax)
        
        self.checkflux(plot=plot,out=out,fluxmin=fluxmin,fluxmax=fluxmax,fitinfo=fitinfo)
        if plot:
            plt.show()

    def calibiodet(self,**kwargs):
        self.caliboption = kwargs.get("caliboption","optics")
        flux = self.idet.cpstoflux(self.energy,self.idetcps)
        #flux = self.fluxest # this is the flux calculated in spec
        
        return self.iodet.calibrate(self.iodetcps,flux,self.energy,**kwargs)
        
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
        self.defaulttime = units.Quantity(time,"s")

    def xrfnormop(self,energy,time=None,ref=None,referencetime=None):
        """
        Returns:
            op(linop): raw diode conversion operator
            Fref(num): flux in photons/s to which the data is normalized after data/op(diode)
            tref(num): time in s to which the data is normalized after data/op(diode)
            t(num):    time in s of raw data
        """
        if time is None:
            time = self.defaulttime
        if ref is None:
            ref = self.reference
        if referencetime is None:
            referencetime = self.referencetime
            
        ret = self.iodet.xrfnormop(energy,time,ref,referencetime=referencetime)
        return ret+(units.magnitude(time,"s"),)
    
    def fluxtocps(self,energy,flux):
        return self.iodet.fluxtocps(energy,flux).magnitude
  
    def cpstoflux(self,energy,cps):
        return self.iodet.cpstoflux(energy,cps).magnitude
    
    def I0op(self,energy,time=None):
        if time is None:
            time = self.defaulttime
        op = self.iodet.fluxop(energy,time)
        return op,time

    def Itop(self,energy,time=None):
        if time is None:
            time = self.defaulttime
        op = self.idet.fluxop(energy,time)
        return op,time
        
    def plot_response(self,energy=None,flux=1e9,diode="iodet",current=False):
        if energy is None:
            energy = np.linspace(1,15,50)
            
        if flux is None:
            func = lambda o: getattr(o, 'plot_spectral_responsivity')(energy)
        else:
            func = lambda o: getattr(o, 'plot_response')(energy,flux,current=current)

        if diode=="iodet":
            func(self.iodet)
        elif diode=="idet":
            func(self.idet)
        elif diode=="ptb":
            func(self.idet.absdiode)
        plt.show()
        
        # Check idet response:
        #fluxmonitor.plot_response(diode='idet',current=True)
        #fluxmonitor.idet.spectral_responsivity = super(spectrocrunch.detectors.diode.CalibratedPNdiode,fluxmonitor.idet).spectral_responsivity
        #fluxmonitor.plot_response(diode='idet',current=True)
        #import matplotlib.pyplot as plt
        #plt.show()
        
    def calibrate(self,params):
        # required
        base = params["base"]
        sample = params["sample"]
        dataset = params["dataset"]
        specnr = params["specnr"]
        dark = params["dark"]
        
        # optional
        gainiodet = params.get("gainiodet",None)
        gainidet = params.get("gainidet",None)
        fixdark = params.get("fixdark",False)
        nofit = params.get("nofit",False)
        plot = params.get("plot",False)
        fluxmin = params.get("fluxmin",0)
        fluxmax = params.get("fluxmax",np.inf)
        
        # map info from spec
        sampledataset = "{}_{}".format(sample,dataset)
        specfile = os.path.join(base,sample,sampledataset,"{}.dat".format(sampledataset))
        
        if dark:
            # loopscan
            nofit |= fixdark
            self.dark(specfile=specfile,specnr=specnr,gainiodet=gainiodet,gainidet=gainidet,nofit=nofit,plot=plot)
        else:
            # attz scan
            self.calib(specfile=specfile,specnr=specnr,gainiodet=gainiodet,gainidet=gainidet,\
                    nofit=nofit,fixdark=fixdark,fluxmin=fluxmin,fluxmax=fluxmax,plot=plot)

    def batchcalibrate(self,params_fixed,params_var):
        for k in params_var:
            params = dict(params_fixed)
            params.update(k)
            self.calibrate(params)
        
        return
        
        #print(monitor)
        #print(monitor.iodet)
        #exit()

        #energy = np.linspace(1,15,50)
        #flux = 1e9
        #plt.plot(energy,np.asarray([monitor.iodet._detection_rates(np.asarray([en]),np.asarray([1]))[1].sum() for en in energy]))
        #plt.plot(energy,[monitor.iodet._chargepersamplephoton(en) for en in energy])
        #self.plot_response()
        
    
        
        

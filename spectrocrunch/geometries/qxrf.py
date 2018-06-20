# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
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

from . import xrf as xrfgeometries
from ..detectors import diode as diodes
from ..detectors import xrf as xrfdetectors
from ..sources import xray as xraysources
from ..optics import xray as xrayoptics
from ..instruments import configuration
from ..common import units
from ..common import instance
from ..common.classfactory import with_metaclass
from ..io import spec
from .. import ureg

import numpy as np
import matplotlib.pyplot as plt
from pint import errors as pinterrors

logger = logging.getLogger(__name__)


class QXRFGeometry(with_metaclass(object)):
    """Quantitative XRF geometry with I0 and It diodes
    """
    
    def __init__(self,instrument=None,diodeI0=None,diodeIt=None,xrfdetector=None,\
                xrfgeometry=None,optics=None,source="synchrotron",simplecalibration=True):
        
        self.instrument = configuration.factory(instrument)

        self._diodeI0 = None
        self._diodeIt = None
        self._xrfgeometry = None
        
        self.simplecalibration = simplecalibration
        self.optics = optics
        self.source = source
        self.xrfdetector = xrfdetector
        self.xrfgeometry = xrfgeometry
        self.diodeI0 = diodeI0
        self.diodeIt = diodeIt

        self.reference = units.Quantity(1e10,"hertz")
        self.defaultexpotime = units.Quantity(0.1,"s")
        self.referencetime = None
        self._data = {}

    def __str__(self):
        if self.diodeI0 is None:
            diodeI0 = None
        else:
            diodeI0 = self.diodeI0.geometry
            if hasattr(diodeI0,"geometry"):
                diodeI0 = diodeI0.geometry
        
        if self.diodeIt is None:
            diodeIt = None
        else:
            diodeIt = self.diodeIt
            if hasattr(diodeIt,"geometry"):
                diodeIt = diodeIt.geometry
        
        referencetime = self.referencetime
        if referencetime is not None:
            referencetime = "{:~}".format(referencetime)

        fmt = "DiodeI0:\n========\n{}\n"\
               "DiodeIt:\n========\n{}\n"\
               "XRF geometry:\n===========\n{}\n"\
               "Flux monitor:\n===========\n"\
               " Default flux reference:\n {:~e}\n"\
               " Reference exposure time:\n {}\n"\
               " Default exposure time:\n {:~}\n"
               
        return fmt.format(diodeI0,diodeIt,self.xrfgeometry,\
                   self.reference,referencetime,self.defaultexpotime)
    
    @property
    def simplecalibration(self):
        return self._simplecalibration

    @simplecalibration.setter
    def simplecalibration(self,value):
        self._simplecalibration = value
        self._update_devices()
        
    @property
    def diodeI0(self):
        return self._diodeI0

    @diodeI0.setter
    def diodeI0(self,device):
        self._diodeI0 = self._generate_device(device,diodes.factory)
        self._update_devices()
    
    @property
    def diodeIt(self):
        return self._diodeIt

    @diodeIt.setter
    def diodeIt(self,device):
        self._diodeIt = self._generate_device(device,diodes.factory)
        self._update_devices()
    
    @property
    def source(self):
        return self._source
        
    @source.setter
    def source(self,device):
        self._source = self._generate_device(device,xraysources.factory)
        self._update_devices()
    
    @property
    def xrfdetector(self):
        return self._xrfdetector
        
    @xrfdetector.setter
    def xrfdetector(self,device):
        self._xrfdetector = self._generate_device(device,xrfdetectors.factory)
        self._update_devices()
    
    @property
    def optics(self):
        return self._optics
        
    @optics.setter
    def optics(self,device):
        self._optics = self._generate_device(device,xrayoptics.factory)
        self._update_devices()
        
    @property
    def xrfgeometry(self):
        return self._xrfgeometry

    @xrfgeometry.setter
    def xrfgeometry(self,device):
        self._xrfgeometry = self._generate_device(device,xrfgeometries.factory)
        self._update_devices()

    def _generate_device(self,device,factory):
        if instance.isstring(device):
            return factory(device)
        else:
            return device
            
    def _update_devices(self):
        if self.diodeI0 is not None:
            self.diodeI0.optics = self.optics
            if self.diodeI0.secondarytarget is not None:
                self.diodeI0.secondarytarget.geometry.source = self.source
            self.diodeI0.simplecalibration = self.simplecalibration
        if self.xrfgeometry is not None:
            self.xrfgeometry.detector = self.xrfdetector
            self.xrfgeometry.source = self.source

    def setxrfposition(self,value):
        self.xrfgeometry.detectorposition = value
    
    def getxrfdistance(self):
        return self.xrfgeometry.distance.to("cm").magnitude

    def setreferenceflux(self,flux):
        self.reference = units.Quantity(flux,"hertz")
    
    def setreferencetime(self,expotime):
        if expotime is None:
            self.referencetime = None # take the time of the scan
        else:
            self.referencetime = units.Quantity(expotime,"s")
        
    def setreferencecounts(self,cts):
        self.reference = units.Quantity(cts,"dimensionless")
    
    def setdefaulttime(self,expotime):
        self.defaultexpotime = units.Quantity(expotime,"s")

    def fluxtocps(self,energy,flux,weights=None):
        return self.diodeI0.fluxtocps(energy,flux,weights=weights).magnitude
  
    def responsetoflux(self,energy,response,weights=None):
        return self.diodeI0.responsetoflux(energy,response,weights=weights).magnitude
    
    def xrfnormop(self,energy,expotime=None,reference=None,referencetime=None,weights=None):
        """
        Args:
            energy(num|array): source lines (keV)
            expotime(Optional(num)): original exposure time (sec)
            reference(Optional(num)): iodet (counts) or flux (photons/sec) to which the data should be normalized
            referencetime(Optional(num)): time to which the data should be normalized
            weights(Optional(num|array)): source line weights
            
        Returns:
            op(linop): raw diode conversion operator
            Fref(num): flux in photons/s to which the data is normalized after data/op(diode)
            tref(num): time in s to which the data is normalized after data/op(diode)
            t(num):    time in s of raw data
        """
        if expotime is None:
            expotime = self.defaultexpotime
        if reference is None:
            reference = self.reference
        if referencetime is None:
            referencetime = self.referencetime
            
        ret = self.diodeI0.xrfnormop(energy,expotime,reference,referencetime=referencetime,weights=weights)
        return ret+(units.umagnitude(expotime,"s"),)
        
    def I0op(self,energy,expotime=None,weights=None):
        if expotime is None:
            expotime = self.defaultexpotime
        op = self.diodeI0.fluxop(energy,expotime)
        return op,expotime

    def Itop(self,energy,expotime=None,weights=None):
        if expotime is None:
            expotime = self.defaultexpotime
        op = self.diodeIt.fluxop(energy,expotime,weights=weights)
        return op,expotime

    def batchcalibrate(self,params_fixed,params_var):
        for k in params_var:
            params = dict(params_fixed)
            params.update(k)
            self.calibrate(**params)
            
    def calibrate(self,**paramsin):
        params = dict(paramsin)
        
        # Get data
        base = params.pop("base",None)
        sample = params.pop("sample",None)
        dataset = params.pop("dataset",None)
        specnr = params.pop("specnr",None)

        if specnr is None:
            # Data in the params
            data,staticdata = self._parse_default(params)
        else:
            sampledataset = "{}_{}".format(sample,dataset)
            specfile = os.path.join(base,sample,sampledataset,"{}.dat".format(sampledataset))
            data,staticdata = self._parse_spec(specfile,specnr)

        # 
        resetdevices = params.pop("resetdevices",False)
        params2 = self._validate_data(data,staticdata,resetdevices=resetdevices)
        
        # Optional parameters
        gaindiodeI0 = params.pop("gaindiodeI0",None)
        gaindiodeIt = params.pop("gaindiodeIt",None)
        nofit = params.pop("nofit",False)
        plot = params.pop("plot",False)
        
        # Set diode gains
        self._set_diode_gain(gaindiodeI0,params2.get("gaindiodeI0",None),"diodeI0")
        self._set_diode_gain(gaindiodeIt,params2.get("gaindiodeIt",None),"diodeIt")

        # Calibrate
        dark = params.pop("dark") 
        if dark:
            fixdark = params.pop("fixdark",False)
            params.pop("fluxmin",None)
            params.pop("fluxmax",None)
            if not(nofit or fixdark):
                self.diodeI0.darkfromcps(data["I0"],**params)
                self.diodeIt.darkfromcps(data["It"],**params)
            if plot:
                self._show_dark_calib(data)
        else:
            if nofit:
                fitinfo = {}
            else:
                energy = data["energy"].to("keV").magnitude
                flux = self.diodeIt.responsetoflux(energy,data["It"])
                fitinfo = self.diodeI0.calibrate(data["I0"],flux,energy,**params)
            if plot:
                fluxmin = params.get("fluxmin",0)
                fluxmax = params.get("fluxmax",np.inf)
                self._show_flux_calib(data,fluxmin=fluxmin,fluxmax=fluxmax,fitinfo=fitinfo)
        if plot:
            plt.show()

    def _set_diode_gain(self,gainasked,gaindata,attr):
        deviceinstance = getattr(self,attr)
        if gainasked is None:
            if gaindata is None:
                raise RuntimeError("Gain of diode {} must be specified".format(deviceinstance.__class__.__name__))
            deviceinstance.gain = gaindata
        else:
            if gaindata is not None:
                gainasked = units.quantity_like(gainasked,gaindata)
                if gainasked!=gaindata:
                    raise RuntimeError("Data was collected with diode {} gain {}, but {} was expected".format(gaindata,deviceinstance.__class__.__name__,gainasked))
            deviceinstance.gain = gainasked

    def _check_device(self,attr,clsfactory,clsnameused,resetdevice=False):
        if clsnameused is None:
            return # usage is not specified
            
        deviceinstance = getattr(self,attr)
        if clsnameused:
            # Device was used
            if deviceinstance is None:
                # Device was not expected to be used
                if resetdevice:
                    setattr(self,attr,clsnameused)
                else:
                    raise RuntimeError("Device {} was not expected to be used.".format(clsnameused))
            else:
                # Device was expected to be used
                deviceclass = clsfactory(clsnameused)
                if not isinstance(deviceinstance,deviceclass):
                    raise RuntimeError("Data was collected with diode {}, while diode {} was expected.".format(deviceclass.__name__,deviceinstance.__class__.__name__))
        else:
            # Device was not used
            if deviceinstance is not None:
                # Device was expected to be used
                if resetdevice:
                    setattr(self,attr,clsnameused)
                else:
                    raise RuntimeError("Data was not collected with device {}.".format(deviceinstance.__class__.__name__))
    
    def _calibrate_fields(self):
        return ["I0","It","flux0","fluxt","time","energy"]
    
    def _parse_spec(self,specfile,specnr):
        fspec = spec.spec(specfile)
        
        # Motors
        motornames = self.instrument.motornames()
        specmotornames = self.instrument.specmotornames(motornames)
        motorvalues = fspec.getmotorvalues(specnr,specmotornames)
        staticdata = dict(zip(motornames,motorvalues))

        # Counters
        fields = self._calibrate_fields()
        labels = [self.instrument.speccounternames.get(k,"") for k in fields]
        available = fspec.haslabels(specnr,labels)
        labels = [label for label,v in zip(labels,available) if v]
        ufields = [name for name,v in zip(fields,available) if v]
        data = fspec.getdata2(specnr,labels).T
        data = {k:units.Quantity(v,self.instrument.units[k]) for k,v in zip(ufields,data)}
        
        return data,staticdata
    
    def _parse_default(self,params):
        fields = self._calibrate_fields()
        
        data = {}
        for k in fields:
            if k in params:
                data[k] = units.Quantity(params.pop(k),self.instrument.units[k])
                
        return data,params.pop("motors",{})
        
    def _validate_data(self,data,staticdata,resetdevices=False):
        
        fields = self._calibrate_fields()
        
        # Add units where missing and fill with static data if needed
        for k in fields:
            if k in data:
                data[k] = units.Quantity(data[k],self.instrument.units[k])
            elif k in staticdata:
                data[k] = units.Quantity(staticdata[k],self.instrument.units[k])
        
        # Fill in missing exposure time
        if "time" not in data:
            data["time"] = self.defaultexpotime
        data["time"] = data["time"].to("s")
        
        # Check units of energy
        if "energy" in data:
            data["energy"] = data["energy"].to("keV")
        
        # Make sure flux is in hertz
        for k in ["flux0","fluxt"]:
            if k in data:
                if data[k].units == ureg.dimensionless:
                    data[k] = data[k]/data["time"]
                data[k] = data[k].to("Hz")
        
        # Make sure diode response is in Hz(default) or A
        for k in ["I0","It"]:
            if k in data:
                if data[k].units == ureg.dimensionless:
                    data[k] = data[k]/data["time"]
                data[k] = units.unitsto(data[k],["Hz","A"])

        # Verify, create or replace devices
        if "I0" in data:
            try:
                diodeI0name = self.instrument.diodeI0(staticdata)
            except KeyError:
                diodeI0name = None
        else:
            diodeI0name = ""
        
        if "It" in data:
            try:
                diodeItname = self.instrument.diodeIt(staticdata)
            except KeyError:
                diodeItname = None
        else:
            diodeItname0name = ""
        
        try:
            opticsname = self.instrument.optics(staticdata)
        except KeyError:
            opticsname = None

        self._check_device("optics",xrayoptics.clsfactory,opticsname,resetdevice=resetdevices)
        self._check_device("diodeI0",diodes.clsfactory,diodeI0name,resetdevice=resetdevices)
        self._check_device("diodeIt",diodes.clsfactory,diodeItname,resetdevice=resetdevices)
        
        # Set gains
        ret = {}
        if "I0" in data and "flux0" in data and "energy" in data:
            ret["gaindiodeI0"] = self.diodeI0.gainfromresponse(data["energy"].to("keV").magnitude,data["I0"],data["flux0"])

        gainIt = "It" in data and "fluxt" in data
        if "It" in data and "fluxt" in data and "energy" in data:
            ret["gaindiodeIt"] = self.diodeIt.gainfromresponse(data["energy"].to("keV").magnitude,data["It"],data["fluxt"])

        return ret

    def _show_dark_calib(self,data):
        f, (ax1, ax2) = plt.subplots(1, 2)
        self._plot_dark_calib(self.diodeIt,data["It"],"diodeIt",ax1)
        color = next(ax1._get_lines.prop_cycler)['color']
        self._plot_dark_calib(self.diodeI0,data["I0"],"diodeI0",ax2,color=color)
        plt.tight_layout()

    def _plot_dark_calib(self,deviceinstance,response,name,ax,color=None):
        ax.plot(response,'o',label='spec',color=color)
        ax.axhline(y=deviceinstance.fluxtocps(5,0).magnitude,label='calc',color=color)
        ax.set_xlabel("points")
        ax.set_ylabel("{} ({:~})".format(name,response.units))
        ax.set_title("{:~.0e}".format(deviceinstance.gain))
        ax.legend(loc="best")

    def _plot_diode_flux(self,attr,energy,response,ax,color=None):
        deviceinstance = getattr(self,attr)
        
        if response.size==1:
            response = units.asqarray([response[0]*0.9,response[0]*1.1])
            
        flux = deviceinstance.responsetoflux(energy,response)
        
        lines = ax.plot(flux,response,label=attr,color=color)
        ax.set_xlabel("flux (ph/s)")
        ax.set_ylabel("{} ({:~})".format(attr,response.units))
        ax.set_title("{:~.0e}".format(deviceinstance.gain))
        ax.legend(loc="best")
        
        return lines[0].get_color()
        
    def _show_flux_calib(self,data,fluxmin=0,fluxmax=np.inf,fitinfo={}):
        f, (ax1, ax2) = plt.subplots(1, 2)
        
        energy = data["energy"].to("keV").magnitude
        
        x = units.asqarray(data["fluxt"])
        It = units.asqarray(data["It"])
        ax1.plot(x,It,'o',label='spec')
        color1 = self._plot_diode_flux('diodeIt',energy,It,ax1)
        color2 = next(ax1._get_lines.prop_cycler)['color']
        
        xoff = 0.3
        yoff = 0.01
        info = self.diodeIt.fluxcpsinfo(energy)
        text = "\n".join(["{} = {}".format(k,v) for k,v in info.items()])
        ax1.text(xoff,yoff,text,transform = ax1.transAxes,verticalalignment="bottom")

        x = self.diodeIt.responsetoflux(energy,It)
        fluxmin = units.Quantity(fluxmin,"hertz").to("hertz")
        fluxmax = units.Quantity(fluxmax,"hertz").to("hertz")
        I0 = units.asqarray(data["I0"])
        indfit = (x<fluxmin) | (x>fluxmax)
        if any(indfit):
            ax2.plot(x[indfit],I0[indfit],'o',color="#cccccc")
        indfit = np.logical_not(indfit)
        if any(indfit):
            ax2.plot(x[indfit],I0[indfit],'o',label='diodeIt',color=color1)
        
        self._plot_diode_flux('diodeI0',energy,I0,ax2,color=color2)

        if fitinfo:
            text = "\n".join(["{} = {}".format(k,v) for k,v in fitinfo.items()])
            ax2.text(xoff,yoff,text,transform = ax2.transAxes,verticalalignment="bottom")
            
        plt.tight_layout()


factory = QXRFGeometry.factory


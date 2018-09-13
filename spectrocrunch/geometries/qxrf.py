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
from ..patch.pint import ureg
from ..math.linop import LinearOperator
from ..math.common import weightedsum

import numpy as np
import matplotlib.pyplot as plt
from pint import errors as pinterrors

logger = logging.getLogger(__name__)


class QXRFGeometry(with_metaclass(object)):
    """Quantitative XRF geometry with I0 and It diodes
    """
    
    def __init__(self,instrument=None,diodeI0=None,diodeIt=None,optics=None,\
                xrfgeometries=None,xrfdetector=None,xrfgeometry=None,source="synchrotron",\
                simplecalibration=True,fluxid="I0",transmissionid="It"):
        
        self.instrument = configuration.factory(instrument)
        
        if xrfgeometries is None:
            xrfgeometries = [(xrfgeometry,xrfdetector)]

        self._diodeI0 = None
        self._diodeIt = None
        self._xrfgeometries = []
        self.simplecalibration = simplecalibration
        self.optics = optics
        self.source = source
        self.xrfgeometries = xrfgeometries
        self.diodeI0 = diodeI0
        self.diodeIt = diodeIt

        self.reference = units.Quantity(1e10,"hertz")
        self.defaultexpotime = units.Quantity(0.1,"s")
        self.referencetime = None
        self._data = {}
        self.fluxid = fluxid
        self.transmissionid = transmissionid

    def __str__(self):
        lst = []
        
        if self.diodeI0 is None:
            diodeI0 = None
        else:
            diodeI0 = self.diodeI0
            if hasattr(diodeI0,"geometry"):
                diodeI0 = diodeI0.geometry
        lst.append("DiodeI0:\n========\n{}".format(diodeI0))
        
        if self.diodeIt is None:
            diodeIt = None
        else:
            diodeIt = self.diodeIt
            if hasattr(diodeIt,"geometry"):
                diodeIt = diodeIt.geometry
        lst.append("DiodeIt:\n========\n{}".format(diodeIt))
        
        for i,geom in enumerate(self.xrfgeometries,1):
            lst.append("XRF geometry ({}):\n===========\n{}".format(i,geom))

        lst.append("Flux monitor:\n===========") 
        try:
            self.reference.to("Hz")
            refname = "Normalize to flux"
        except pinterrors.DimensionalityError:
            refname = "Normalize to diode counts"
        lst.append("{}:\n {:~e}".format(refname,self.reference))
        referencetime = self.referencetime
        if referencetime is None:
            referencetime = "<data>"
        else:
            referencetime = "{:~}".format(referencetime)
        lst.append("Normalize to exposure time:\n {}".format(referencetime))
        lst.append("Default exposure time:\n {:~}",self.defaultexpotime)
        
        return '\n'.join(lst)
    
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
    def optics(self):
        return self._optics
        
    @optics.setter
    def optics(self,device):
        self._optics = self._generate_device(device,xrayoptics.factory)
        self._update_devices()
        
    @property
    def xrfgeometries(self):
        return self._xrfgeometries

    @property
    def xrfgeometry(self):
        if len(self._xrfgeometries)!=1:
            raise RuntimeError('More than one XRF geometry, use "xrfgeometries"')
        return self.xrfgeometries[0]

    @xrfgeometries.setter
    def xrfgeometries(self,detgeoms):
        self._xrfgeometries = []
        if detgeoms:
            for geom in detgeoms:
                if not isinstance(geom,xrfgeometries.XRFGeometry):
                    xrfgeometry,xrfdetector = geom
                    geom = self._generate_device(xrfgeometry,xrfgeometries.factory)
                    if geom:
                        geom.detector = self._generate_device(xrfdetector,xrfdetectors.factory)
                if geom:
                    self._xrfgeometries.append(geom)
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
        if self.xrfgeometries and self.source is not None:
            for geom in self.xrfgeometries:
                geom.source = self.source

    def setxrfposition(self,value):
        if instance.isnumber(value):
            value = [value]
        for geom,pos in zip(self.xrfgeometries,value):
            geom.detectorposition = value
    
    def getxrfdistance(self):
        return [geom.distance.to("cm").magnitude for geom in self.xrfgeometries]

    def getxrfdetectorposition(self):
        return [geom.detectorposition.magnitude for geom in self.xrfgeometries]

    def getxrfactivearea(self):
        return [geom.detector.activearea.to("cm**2").magnitude for geom in self.xrfgeometries]

    def getxrfanglein(self):
        return [geom.anglein for geom in self.xrfgeometries]

    def getxrfangleout(self):
        return [geom.angleout for geom in self.xrfgeometries]
        
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
  
    def responsetoflux(self,energy,response,weights=None,time=None):
        return self.diodeI0.responsetoflux(energy,response,weights=weights,time=time).magnitude
    
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
        else:
            expotime = units.Quantity(expotime,"s")
            
        if reference is None:
            reference = self.reference
        else:
            reference = units.Quantity(reference,"dimensionless")
            
        if referencetime is None:
            referencetime = self.referencetime
        else:
            referencetime = units.Quantity(referencetime,"s")
            
        ret = self.diodeI0.xrfnormop(energy,expotime,reference,referencetime=referencetime,weights=weights)
        return ret+(units.umagnitude(expotime,"s"),)
    
    def quantinfo(self,*args,**kwargs):
        """
        Args:
            see xrfnormop
        
        Returns:
            op(linop): raw diode conversion operator
            info(list(dict)): pymca quantification info (flux and geometry)
        """
        xrfnormop,flux,time,expotime = self.xrfnormop(*args,**kwargs)

        quants = []
        for geom in self.xrfgeometries:
            quant = {'flux':flux,'time':time}
            quant["activearea"] = geom.detector.activearea
            quant["anglein"] = geom.anglein
            quant["angleout"] = geom.angleout
            quant["distance"] = self.getxrfdistance()
            quants.append(quant)
            
        return xrfnormop,quants
    
    def _beamfilter_transmission(self,energy,weights=None):
        Tbeamfilters = [weightedsum(geom.detector.filter_transmission(energy,source=True),weights=weights) for geom in self.xrfgeometries]
        if len(set(Tbeamfilters))>1:
            beamfilters = '\n'.join([' Geometry {}: {}'.format(i,geom.beamfilters()) for i,geom in enumerate(self.xrfgeometries,1)])
            raise RuntimeError('Geometries should have the same beam filters:\n{}'.format(beamfilters))
        else:
            Tbeamfilters = Tbeamfilters[0]
        return Tbeamfilters
        
    def I0op(self,energy,expotime=None,weights=None,removebeamfilters=False):
        """Calculate the flux before the sample from the I0 diode response
        
        Args:
            energy(num|array): keV
            expotime(Optional(num)): sec
            weights(Optional(num|array)): source lines ratio's
            removebeamfilters(Optional(bool)): remove beam filter transmission
        """
        if expotime is None:
            expotime = self.defaultexpotime
        op = self.diodeI0.fluxop(energy,0,time=expotime,weights=weights)
        
        if removebeamfilters:
            Tbeamfilters = self._beamfilter_transmission(energy,weights=weights)
            if Tbeamfilters!=1:
                op = LinearOperator(1./Tbeamfilters,0)*op

        return op,expotime

    def Itop(self,energy,expotime=None,weights=None,removebeamfilters=False):
        """Calculate the flux after the sample from the It diode response
        
        Args:
            energy(num|array): keV
            expotime(Optional(num)): sec
            weights(Optional(num|array)): source lines ratio's
            removebeamfilters(Optional(bool)): remove beam filter transmission
        """
        if expotime is None:
            expotime = self.defaultexpotime
        op = self.diodeIt.fluxop(energy,0,time=expotime,weights=weights)
        
        if removebeamfilters:
            # Filters after the sample (like atmosphere) are already taken into account
            Tbeamfilters = self._beamfilter_transmission(energy,weights=weights)
            if Tbeamfilters!=1:
                op = LinearOperator(1./Tbeamfilters,0)*op
                       
        return op,expotime

    def batchcalibrate(self,params_fixed,params_var):
        for k in params_var:
            params = dict(params_fixed)
            params.update(k)
            self.calibrate(**params)
    
    def _calibrate_prepare(self,params):
        # Get data
        base = params.pop("base",None)
        sample = params.pop("sample",None)
        dataset = params.pop("dataset",None)
        specnr = params.pop("specnr",None)
        
        if specnr:
            if sample:
                sampledataset = "{}_{}".format(sample,dataset)
                specfile = os.path.join(base,sample,sampledataset,"{}.dat".format(sampledataset))
            else:
                specfile = os.path.join(base,dataset)

            # Data in spec file
            data,staticdata = self._parse_spec(specfile,specnr)
        else:
            # Data in params
            data,staticdata = self._parse_default(params)

        # Prepare devices
        resetdevices = params.pop("resetdevices",False)
        params2 = self._validate_data(data,staticdata,resetdevices=resetdevices,dark=params["dark"])
        gaindiodeI0 = params.pop("gaindiodeI0",None)
        gaindiodeIt = params.pop("gaindiodeIt",None)
        self._set_diode_gain(gaindiodeI0,params2.get("gaindiodeI0",None),"diodeI0")
        self._set_diode_gain(gaindiodeIt,params2.get("gaindiodeIt",None),"diodeIt")

        return data

    def _calibrate_dark(self,params,data):
        docalibrate = params.pop("calibrate",True)
        fixdark = params.pop("fixdark",False)
        params.pop("fluxmin",None)
        params.pop("fluxmax",None)
        if docalibrate and not fixdark:
           for k,attr in zip(self._calibrate_fields_id(),["diodeI0","diodeIt"]):
                deviceinstance = getattr(self,attr)
                kcurrent = "{}_current".format(k)
                kcps = "{}_cps".format(k)
                if kcurrent in data and kcps in data and False:
                    deviceinstance.calibratedark(data[kcurrent])
                    deviceinstance.calibrateF0(data[kcps])
                elif kcps in data:
                    deviceinstance.calibratedark(data[kcps])

    def _calibrate_gain(self,params,data):
        fitinfo = {}

        docalibrate = params.pop("calibrate",True)
        if docalibrate:
            energy = data["energy"].to("keV").magnitude
            for name in ["current","cps"]:
                I0,It = self._calibrate_fields_name(name)
                if It in data and I0 in data:
                    flux = self.diodeIt.responsetoflux(energy,data[It])
                    fitinfo = self.diodeI0.calibrate(data[I0],flux,energy,**params)
                    break

        return fitinfo

    def calibrate(self,**paramsin):
        # Get data and prepare devices
        params = dict(paramsin)
        data = self._calibrate_prepare(params)

        # Calibrate
        dark = params.pop("dark",False)
        plot = params.pop("plot",False)
        oscillator = params.pop("oscillator",False)
        if oscillator:
            if plot:
                self._show_oscillator_calib(data)
        else:
            if dark:
                self._calibrate_dark(params,data)
                if plot:
                    self._show_dark_calib(data)
            else:
                fitinfo = self._calibrate_gain(params,data)
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
                if not np.isclose(gainasked,gaindata):
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
        elif clsnameused is None:
            # No idea whether device was used or not
            pass
        else:
            # Device was not used
            if deviceinstance is not None:
                # Device was expected to be used
                if resetdevice:
                    setattr(self,attr,clsnameused)
                else:
                    raise RuntimeError("Data was not collected with device {}.".format(deviceinstance.__class__.__name__))
    
    def _calibrate_fields_id(self):
        return [self.fluxid,self.transmissionid]

    def _calibrate_fields_name(self,name):
        return ["{}_{}".format(k,name) for k in self._calibrate_fields_id()]

    def _calibrate_fields(self):
        return self._calibrate_fields_name("counts")+\
               self._calibrate_fields_name("cps")+\
               self._calibrate_fields_name("current")+\
               self._calibrate_fields_name("photons")+\
               self._calibrate_fields_name("flux")+\
               ["time","energy"]
    
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
        ascounter = fspec.haslabels(specnr,labels)
        asmotor = fspec.hasmotors(specnr,labels)

        labels1 = [label for label,v in zip(labels,ascounter) if v]
        ufields1 = [name for name,v in zip(fields,ascounter) if v]
        data = fspec.getdata2(specnr,labels1).T
        data = {k:units.Quantity(v,self.instrument.units[k]) for k,v in zip(ufields1,data)}
        
        labels2 = [label for label,v in zip(labels,asmotor) if v]
        ufields2 = [name for name,v in zip(fields,asmotor) if v]
        motorvalues = fspec.getmotorvalues(specnr,labels2)
        
        data.update({k:units.Quantity(v,self.instrument.units[k]) for k,v in zip(ufields2,motorvalues)})
        return data,staticdata
    
    def _parse_default(self,params):
        fields = self._calibrate_fields()
        
        data = {}
        for k in fields:
            if k in params:
                data[k] = units.Quantity(params.pop(k),self.instrument.units[k])
                
        return data,params.pop("motors",{})
        
    def _validate_data(self,data,staticdata,resetdevices=False,dark=False):

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
        
        # Diode response: ensure units
        for k in self._calibrate_fields_id():
            # counts are not used anywhere, just needed to calculate cps
            kcounts = "{}_counts".format(k)
            kphotons = "{}_photons".format(k)
            kcurrent = "{}_current".format(k)
            kcps = "{}_cps".format(k)
            kflux = "{}_flux".format(k)

            if kcounts in data:
                data[kcounts] = data[kcounts].to("dimensionless")
                if kcps not in data:
                    data[kcps] = data[kcounts].to("dimensionless")/data["time"]

            if kcps in data:
                data[kcps] = data[kcps].to("Hz")

            if kphotons in data:
                data[kphotons] = data[kphotons].to("dimensionless")
                if kflux not in data:
                    data[kflux] = data[kphotons].to("dimensionless")/data["time"]

            if kflux in data:
                data[kflux] = data[kflux].to("Hz")

            if kcurrent in data:
                data[kcurrent] = data[kcurrent].to("A")

        # Verify, create or replace devices
        I0,It = self._calibrate_fields_name("cps")
        diodeI0name = self.instrument.diodeI0(staticdata)
        diodeItname = self.instrument.diodeIt(staticdata)
        opticsname = self.instrument.optics(staticdata)
        self._check_device("optics",xrayoptics.clsfactory,opticsname,resetdevice=resetdevices)
        self._check_device("diodeI0",diodes.clsfactory,diodeI0name,resetdevice=resetdevices)
        self._check_device("diodeIt",diodes.clsfactory,diodeItname,resetdevice=resetdevices)
        
        # Diode gains from response
        ret = {}
        if "energy" in data and not dark:
            energy = data["energy"].to("keV").magnitude
            for k,attr in zip(self._calibrate_fields_id(),["diodeI0","diodeIt"]):
                kcps = "{}_cps".format(k)
                if kcps in data:
                    for name in ["current","flux"]:
                        j = "{}_{}".format(k,name)
                        if j in data:
                            ret["gain"+attr] = getattr(self,attr).gainfromresponse(data[kcps],data[j],energy=energy)
                            break

        return ret

    def _show_dark_calib(self,data):
        f, axs = plt.subplots(1, 2)
        color = None
        for k,attr,ax in zip(self._calibrate_fields_name("cps"),["diodeI0","diodeIt"],axs):
            self._plot_dark_calib(attr,data[k],attr,ax,color=color)
            color = next(ax._get_lines.prop_cycler)['color']
        plt.tight_layout()

    def _plot_dark_calib(self,deviceattr,response,name,ax,color=None):
        deviceinstance = getattr(self,deviceattr)
        ax.plot(response,'o',label='spec',color=color)
        ax.axhline(y=deviceinstance.fluxtocps(5,0).magnitude,label='calc',color=color)
        ax.set_xlabel("points")
        ax.set_ylabel("{} ({:~})".format(name,response.units))
        ax.set_title("{:~.1e}".format(deviceinstance.gain))
        ax.legend(loc="best")

    def _plot_diode_flux(self,attr,energy,response,ax,color=None):
        deviceinstance = getattr(self,attr)
        
        if response.size==1:
            response = units.asqarray([response[0]*0.9,response[0]*1.1])
            
        flux = deviceinstance.responsetoflux(energy,response)
        
        lines = ax.plot(flux,response,label=attr,color=color)
        ax.set_xlabel("flux@sample (ph/s)")
        ax.set_ylabel("{} ({:~})".format(attr,response.units))
        ax.set_title("{:~.1e}".format(deviceinstance.gain))
        ax.legend(loc="best")
        
        return lines[0].get_color()
        
    def _show_flux_calib(self,data,fluxmin=0,fluxmax=np.inf,fitinfo={}):
        f, (ax1, ax2) = plt.subplots(1, 2)
        
        I0,It = self._calibrate_fields_name("cps")
        _,fluxt = self._calibrate_fields_name("flux")

        energy = data["energy"].to("keV").magnitude
        It = units.asqarray(data[It])
        if fluxt in data:
            flux = units.asqarray(data[fluxt])
            ax1.plot(flux,It,'o',label='spec')
        color1 = self._plot_diode_flux('diodeIt',energy,It,ax1)
        color2 = next(ax1._get_lines.prop_cycler)['color']
        
        xoff = 0.3
        yoff = 0.01
        info = self.diodeIt.fluxcpsinfo(energy)
        text = "\n".join(["{} = {}".format(k,v) for k,v in info.items()])
        ax1.text(xoff,yoff,text,transform = ax1.transAxes,verticalalignment="bottom")

        flux = self.diodeIt.responsetoflux(energy,It)
        fluxmin = units.Quantity(fluxmin,"hertz").to("hertz")
        fluxmax = units.Quantity(fluxmax,"hertz").to("hertz")
        I0 = units.asqarray(data[I0])
        indfit = (flux<fluxmin) | (flux>fluxmax)
        if any(indfit):
            ax2.plot(flux[indfit],I0[indfit],'o',color="#cccccc")
        indfit = np.logical_not(indfit)
        if any(indfit):
            ax2.plot(flux[indfit],I0[indfit],'o',label='diodeIt',color=color1)
        
        self._plot_diode_flux('diodeI0',energy,I0,ax2,color=color2)

        if fitinfo:
            text = "\n".join(["{} = {}".format(k,v) for k,v in fitinfo.items()])
            ax2.text(xoff,yoff,text,transform = ax2.transAxes,verticalalignment="bottom")
            
        plt.tight_layout()

    def _plot_oscillator(self,attr,current,ax,color=None):
        deviceinstance = getattr(self,attr)

        op = deviceinstance.op_currenttocps()
        cps = op(current).to("Hz")
        current = current.to("pA")

        lines = ax.plot(current,cps,label=attr,color=color)
        ax.set_xlabel("{} ({:~})".format(attr,current.units))
        ax.set_ylabel("{} ({:~})".format(attr,cps.units))
        ax.set_title("{:~.1e}, {:~.1e}".format(op.m.to("Hz/pA"),op.b.to("Hz")))
        ax.legend(loc="best")
        
        return lines[0].get_color()

    def _show_oscillator_calib(self,data):
        f, (ax1, ax2) = plt.subplots(1, 2)
        
        I0,It = self._calibrate_fields_name("cps")
        cur0,curt = self._calibrate_fields_name("current")
        I0,It = data[I0],data[It]
        cur0,curt = data[cur0],data[curt]

        lines = ax1.plot(curt.to("pA"),It.to("Hz"),'o',label='spec')
        color1 = lines[0].get_color()
        color2 = next(ax1._get_lines.prop_cycler)['color']
        self._plot_oscillator('diodeIt',curt,ax1,color=color1)

        ax2.plot(cur0.to("pA"),I0.to("Hz"),'o',label='spec',color=color2)
        self._plot_oscillator('diodeI0',cur0,ax2,color=color2)

    def list_detectors(self):
        return xrfdetectors.aliases.keys()

    def list_sources(self):
        return xraysources.aliases.keys()
        
    def list_optics(self):
        return xrayoptics.aliases.keys()
        
    def list_geometries(self):
        return xrfgeometries.aliases.keys()
        
    def list_diodes(self):
        return diodes.aliases.keys()
    
    def list_instruments(self):
        return configuration.aliases.keys()

    def list_devices(self):
        print("instrument = {}".format(self.list_instruments()))
        print("xrfgeometry = {}".format(self.list_geometries()))
        print("xrfdetector = {}".format(self.list_detectors()))
        print("diodeI0,diodeIt = {}".format(self.list_diodes()))
        print("source = {}".format(self.list_sources()))
        print("optics = {}".format(self.list_optics()))
        
factory = QXRFGeometry.factory


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

from ..materials import compoundfromname
from . import base
from ..materials import element
from ..common.classfactory import with_metaclass
from ..common import constants
from ..common import instance

import numpy as np
import fisx

class Detector(with_metaclass(base.AttenuatorBase)):

    def __init__(self,mcazero=None,mcagain=None,mcanoise=None,mcafano=None,ehole=None,\
                      activearea=None,geometry=None,**kwargs):
        self.mcazero = mcazero # keV
        self.mcagain = mcagain # keV/bin
        self.activearea = activearea # cm^2
        self.geometry = geometry
        self.mcanoise = mcanoise # FWHM in keV
        self.mcafano = mcafano
        self.ehole = ehole
        super(Detector,self).__init__(**kwargs)
        
    def __getattr__(self, attr):
        # when Detector does not have this attribute or property, try passing it on to the geometry
        return getattr(self.geometry,attr)

    def __str__(self):
        return "Solid angle = 4*pi*{} srad\nActive area = {} cm^2\n{}\n{}".format(self.solidangle,self.activearea,str(self.geometry),super(Detector,self).__str__())
    
    def linewidth(self,energy):
        """FWHM in keV
        """
        return np.sqrt(self.mcanoise**2 + self.mcafano*self.ehole.to("keV").magnitude*energy*8*np.log(2))
    
    @property
    def solidangle(self):
        distance = self.distance
        r2 = self.activearea/np.pi # squared disk radius
        return 2.*np.pi*(1.-(distance/np.sqrt(r2+distance**2.)))

    @solidangle.setter
    def solidangle(self,value):
        solidanglefrac = value/(4*np.pi)
        if solidanglefrac>=0.5:
            raise ValueError("Solid angle must be < 2.pi")
        r2 = self.activearea/np.pi # squared disk radius
        self.distance = np.sqrt(r2)*(0.5-solidanglefrac)/np.sqrt((1-solidanglefrac)*solidanglefrac)
        
    def addtoconfig(self,config,energy): 
        super(Detector,self).addtoconfig(config)
        
        config["concentrations"]["area"] = self.activearea
        config["concentrations"]["distance"] = self.distance
        
        mcazero = self.mcazero
        mcagain = self.mcagain
        mcanoise = self.mcanoise
        # Compensate fano for difference in Ehole with pymca's 3.85 eV (ClassMcaTheory)
        mcafano = self.mcafano * self.ehole.to("eV").magnitude / 3.85
        config["detector"]["zero"] = mcazero
        config["detector"]["gain"] = mcagain
        config["detector"]["noise"] = mcanoise
        config["detector"]["fano"] = mcafano
        
        emin = 0.9
        emax = energy+0.5
        xmin = int((emin-mcazero)/mcagain)
        xmax = int((emax-mcazero)/mcagain)
        config["fit"]["xmin"] = xmin
        config["fit"]["xmax"] = xmax
        config["fit"]["use_limit"] = 1
 
    def loadfromconfig(self,config):
        super(Detector,self).loadfromconfig(config)
        
        self.activearea = config["concentrations"]["area"]
        self.distance = config["concentrations"]["distance"]
        self.mcazero = config["detector"]["zero"]
        self.mcagain = config["detector"]["gain"]
        self.mcanoise = config["detector"]["noise"]
        # Compensate fano for difference in Ehole with pymca's 3.85 eV (ClassMcaTheory)
        self.mcafano = config["detector"]["fano"] * 3.85 / self.ehole.to("eV").magnitude

    def addtofisx(self,setup,cfg):
        super(Detector,self).addtofisx(setup,cfg)
        
        if "Detector" in self.attenuators:
            detmaterial = self.attenuators["Detector"]
            material = detmaterial["material"]
            thickness = detmaterial["thickness"]
            detector = fisx.Detector(material.name, material.density, thickness)
            detector.setActiveArea(self.activearea)
            detector.setDistance(self.distance)
            #detector.setMaximumNumberOfEscapePeaks(0)
            setup.setDetector(detector)
        
        self.geometry.addtofisx(setup,cfg)
        
    def absorbance(self,energy):
        if "Detector" in self.attenuators:
            det = self.attenuators["Detector"]
            return det["material"].mass_att_coeff(energy)*(det["thickness"]*det["material"].density)
        else:
            return np.full_like(3,np.inf,dtype=float)

    def transmission(self,energy):
        return np.exp(-self.absorbance(energy))

    def attenuation(self,energy):
        return 1-self.transmission(energy)
    
    def efficiency(self,energysource,energydet):
        """Detector efficiency = S/cos(ain)*T(energysource)*T(energydet)*A(energydet)
            S: solid angle detector
            ain: angle of beam with surface normal
            T: transmission by filters (before sample and detector)
            A: attenuation by the detector crystal
            
        Args:
            energysource: n0
            energydet: n1

        Returns:
            array: n0 x n1
        """
        energysource = instance.asarray(energysource)
        energydet = instance.asarray(energydet)
        
        g = self.solidangle/self.cosnormin
        T0 = super(Detector,self).filter_transmission(energysource,source=True)
        T1 = super(Detector,self).filter_transmission(energydet,source=False)
        A = self.attenuation(energydet)
        
        # the cosine term is put here for convenience (comes from integration over sample thickness)
        
        return (g*T0)[:,np.newaxis]*(T1*A)[np.newaxis,:]
        
class Leia(Detector):
    aliases = ["SGX80"]
    
    def __init__(self,**kwargs):
        ultralene = compoundfromname.compoundfromname("ultralene")
        
        attenuators = {}
        attenuators["FoilDetector"] = {"material":ultralene,"thickness":4e-4} #cm
        attenuators["WindowDetector"] = {"material":element.Element('Be'),"thickness":25e-4}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":450e-4}
        
        kwargs["activearea"] = 80e-2 # cm^2
        kwargs["mcazero"] = 0. # keV
        kwargs["mcagain"] = 5e-3 # keV
        kwargs["mcanoise"] = 0.1 # keV
        kwargs["mcafano"] = 0.114
        kwargs["ehole"] = constants.eholepair_si()

        super(Leia,self).__init__(attenuators=attenuators,**kwargs)

class BB8(Detector):
    aliases = ["SGX50"]
    
    def __init__(self,**kwargs):
        ultralene = compoundfromname.compoundfromname("ultralene")
        
        attenuators = {}
        attenuators["FoilDetector"] = {"material":ultralene,"thickness":4e-4} #cm
        attenuators["WindowDetector"] = {"material":element.Element('Be'),"thickness":12.5e-4}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":450e-4}
  
        kwargs["activearea"] = 50e-2 # cm^2
        kwargs["mcazero"] = 0. # keV
        kwargs["mcagain"] = 5e-3 # keV
        kwargs["mcanoise"] = 0.1 # keV
        kwargs["mcafano"] = 0.114
        kwargs["ehole"] = constants.eholepair_si()
        
        super(BB8,self).__init__(attenuators=attenuators,**kwargs)

factory = Detector.factory


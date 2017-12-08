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
from ..materials import pymca
from ..materials import element
from ..common.classfactory import with_metaclass

import numpy as np

class Detector(with_metaclass(pymca.PymcaAttenuators)):

    def __init__(self,mcazero=None,mcagain=None,activearea=None,geometry=None,**kwargs):
        self.mcazero = mcazero # keV
        self.mcagain = mcagain # keV/bin
        self.activearea = float(activearea) # cm^2
        self.geometry = geometry
        super(Detector,self).__init__(**kwargs)
        
    def __getattr__(self, attr):
        # when Detector does not have this attribute or property, try passing it on to the geometry
        return getattr(self.geometry,attr)

    def __str__(self):
        return "Solid angle = 4*pi*{} srad\n {}".format(self.solidangle,str(self.geometry))
    
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
        config["detector"]["zero"] = mcazero
        config["detector"]["gain"] = mcagain
        
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
 
class Leia(Detector):
    aliases = ["SGX80"]
    
    def __init__(self,**kwargs):
        ultralene = compoundfromname.compoundfromname("ultralene")
        
        attenuators = {}
        attenuators["FoilDetector"] = {"material":ultralene,"thickness":4e-4}
        attenuators["WindowDetector"] = {"material":element.Element('Be'),"thickness":25e-4}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":450e-4}

        super(Leia,self).__init__(activearea=80e-2,attenuators=attenuators,**kwargs)

class BB8(Detector):
    aliases = ["SGX50"]
    
    def __init__(self,**kwargs):
        ultralene = compoundfromname.compoundfromname("ultralene")
        
        attenuators = {}
        attenuators["FoilDetector"] = {"material":ultralene,"thickness":4e-4}
        attenuators["WindowDetector"] = {"material":element.Element('Be'),"thickness":12.5e-4}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":450e-4}
  
        super(BB8,self).__init__(activearea=50e-2,attenuators=attenuators,**kwargs)

factory = Detector.factory


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

import numpy as np

class Detector(with_metaclass(base.PointSourceCentric)):

    def __init__(self,mcazero=None,mcagain=None,mcanoise=None,mcafano=None,**kwargs):
        self.mcazero = mcazero # keV
        self.mcagain = mcagain # keV/bin
        self.mcanoise = mcanoise # FWHM in keV
        self.mcafano = mcafano
        super(Detector,self).__init__(**kwargs)

    def linewidth(self,energy):
        """FWHM in keV
        """
        return np.sqrt(self.mcanoise**2 + self.mcafano*self.ehole.to("keV").magnitude*energy*8*np.log(2))
        
    def addtopymcaconfig(self,config,energy): 
        super(Detector,self).addtopymcaconfig(config)
        
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
 
    def loadfrompymcaconfig(self,config):
        super(Detector,self).loadfrompymcaconfig(config)
        
        self.activearea = config["concentrations"]["area"]
        self.distance = config["concentrations"]["distance"]
        self.mcazero = config["detector"]["zero"]
        self.mcagain = config["detector"]["gain"]
        self.mcanoise = config["detector"]["noise"]
        # Compensate fano for difference in Ehole with pymca's 3.85 eV (ClassMcaTheory)
        self.mcafano = config["detector"]["fano"] * 3.85 / self.ehole.to("eV").magnitude
        
        
class Leia(Detector):
    aliases = ["SGX80"]
    
    def __init__(self,**kwargs):
        ultralene = compoundfromname.compoundfromname("ultralene")
        
        attenuators = {}
        attenuators["FoilDetector"] = {"material":ultralene,"thickness":4e-4} #cm
        attenuators["WindowDetector"] = {"material":element.Element('Be'),"thickness":25e-4}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":450e-4}
        kwargs["attenuators"] = attenuators
        
        kwargs["activearea"] = 80e-2 # cm^2
        kwargs["mcazero"] = 0. # keV
        kwargs["mcagain"] = 5e-3 # keV
        kwargs["mcanoise"] = 0.1 # keV
        kwargs["mcafano"] = 0.114
        kwargs["ehole"] = constants.eholepair_si()

        super(Leia,self).__init__(**kwargs)


class BB8(Detector):
    aliases = ["SGX50"]
    
    def __init__(self,**kwargs):
        ultralene = compoundfromname.compoundfromname("ultralene")
        
        attenuators = {}
        attenuators["FoilDetector"] = {"material":ultralene,"thickness":4e-4} #cm
        attenuators["WindowDetector"] = {"material":element.Element('Be'),"thickness":12.5e-4}
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":450e-4}
        kwargs["attenuators"] = attenuators
        
        kwargs["activearea"] = 50e-2 # cm^2
        kwargs["mcazero"] = 0. # keV
        kwargs["mcagain"] = 5e-3 # keV
        kwargs["mcanoise"] = 0.1 # keV
        kwargs["mcafano"] = 0.114
        kwargs["ehole"] = constants.eholepair_si()
        
        super(BB8,self).__init__(**kwargs)

factory = Detector.factory


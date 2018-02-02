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

from . import base
from ..materials import compoundfromname
from ..materials import element
from ..common import constants
from ..common.classfactory import with_metaclass

import scipy.special
import numpy as np

class Detector(with_metaclass(base.PointSourceCentric)):

    FWHMlorentz = 0.01 # keV

    def __init__(self,mcazero=None,mcagain=None,mcanoise=None,mcafano=None,\
                    tailwidth=None,tailfraction=None,stepfraction=None,**kwargs):
        self.mcazero = mcazero # keV
        self.mcagain = mcagain # keV/bin
        self.mcanoise = mcanoise # FWHM in keV
        self.mcafano = mcafano
        self.tailwidth = tailwidth
        self.tailfraction = tailfraction
        self.stepfraction = stepfraction
        
        super(Detector,self).__init__(**kwargs)

    def __str__(self):
        return "XRF detector:\n{}".format(super(Detector,self).__str__())

    def linesigma2(self,energy):
        return self.mcanoise**2/(8*np.log(2)) + self.mcafano*self.ehole*1e3*energy

    def lineprofile(self,x,u,voigt=False):
        """
        Args:
            x: energies (keV, nx)
            u: peak energies (keV, nu)
            
        Returns:
            array: nx x nu
        """
        x = instance.asarray(x)[:,np.newaxis]
        u = instance.asarray(u)[np.newaxis,:]

        sigma2 = self.linesigma2(u)
        diff = x-u
        
        a = 2*sigma2
        b = np.sqrt(a)
        c = self.tailwidth*np.sqrt(sigma2)
        detnorm = np.sqrt(np.pi*a)
        
        if voigt:
            detresponse = np.real(scipy.special.wofz((diff + 0.5j*FWHMlorentz)/b))/detnorm
        else:
            detresponse = np.exp(-diff**2/a)/detnorm

        argstep = diff/b
        step = scipy.special.erfc(argstep)/(2*u)
        
        tail = np.exp(1/(2*self.tailwidth**2))/(2*c)*np.exp(diff/c)*scipy.special.erfc(argstep+1/(np.sqrt(2)*self.tailwidth))

        return detrespons + self.tailfraction*tail + self.stepfraction*step
    
    def addtopymca(self,setup,cfg): 
        super(Detector,self).addtopymca(setup,cfg)

        mcazero = self.mcazero
        mcagain = self.mcagain
        mcanoise = self.mcanoise
        # Compensate fano for difference in Ehole with pymca's 3.85 eV (ClassMcaTheory)
        mcafano = self.mcafano * self.ehole / 3.85
        cfg["detector"]["zero"] = mcazero
        cfg["detector"]["gain"] = mcagain
        cfg["detector"]["noise"] = mcanoise
        cfg["detector"]["fano"] = mcafano

        xmin = int((setup.emin-mcazero)/mcagain)
        xmax = int((setup.emax-mcazero)/mcagain)
        cfg["fit"]["xmin"] = xmin
        cfg["fit"]["xmax"] = xmax
        cfg["fit"]["use_limit"] = 1
 
    def loadfrompymca(self,cfg):
        super(Detector,self).loadfrompymca(cfg)
        
        self.activearea = cfg["concentrations"]["area"]
        self.distance = cfg["concentrations"]["distance"]
        self.mcazero = cfg["detector"]["zero"]
        self.mcagain = cfg["detector"]["gain"]
        self.mcanoise = cfg["detector"]["noise"]
        # Compensate fano for difference in Ehole with pymca's 3.85 eV (ClassMcaTheory)
        self.mcafano = cfg["detector"]["fano"] * 3.85 / self.ehole
        

class sn3102(Detector):
    aliases = ["XFlash5100"]
    
    def __init__(self,**kwargs):
        ultralene = compoundfromname.compoundfromname("ultralene")
        moxtek = compoundfromname.compoundfromname("moxtek ap3.3")
        
        attenuators = {}
        attenuators["FoilDetector"] = {"material":ultralene,"thickness":4e-4} # cm
        attenuators["WindowDetector"] = {"material":moxtek,"thickness":380e-4} # AP3.3 Moxtek
        attenuators["Detector"] = {"material":element.Element('Si'),"thickness":450e-4}
        kwargs["attenuators"] = attenuators
        
        kwargs["activearea"] = 80e-2*0.75 # cm^2
        
        kwargs["mcazero"] = 0. # keV
        kwargs["mcagain"] = 5e-3 # keV
        kwargs["mcanoise"] = 0.1 # keV
        kwargs["mcafano"] = 0.114
        kwargs["tailwidth"] = 2.5
        kwargs["tailfraction"] = 0.02 # 2% of the photons end up in the tail
        kwargs["stepfraction"] = 0.03 # 3% of the photons end up in the step
        
        kwargs["ehole"] = constants.eholepair_si()

        super(sn3102,self).__init__(**kwargs)
        
              
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
        kwargs["tailwidth"] = 2.5
        kwargs["tailfraction"] = 0.02 # 2% of the photons end up in the tail
        kwargs["stepfraction"] = 0.03 # 3% of the photons end up in the step
        
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
        kwargs["tailwidth"] = 2.5
        kwargs["tailfraction"] = 0.02 # 2% of the photons end up in the tail
        kwargs["stepfraction"] = 0.03 # 3% of the photons end up in the step
        
        kwargs["ehole"] = constants.eholepair_si()
        
        super(BB8,self).__init__(**kwargs)


class DR40(Detector):
    aliases = ["VITUSH80"]
    
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
        kwargs["tailwidth"] = 2.5
        kwargs["tailfraction"] = 0.02 # 2% of the photons end up in the tail
        kwargs["stepfraction"] = 0.03 # 3% of the photons end up in the step
        
        kwargs["ehole"] = constants.eholepair_si()
        
        super(DR40,self).__init__(**kwargs)
        
factory = Detector.factory


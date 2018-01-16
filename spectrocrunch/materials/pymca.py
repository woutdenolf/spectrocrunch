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

from ..common import units
from . import mixture
from . import compoundfromformula

import fisx
import contextlib
import numpy as np

class PymcaAttenuators(object):

    def __init__(self,attenuators=None):
        self.attenuators = attenuators

    def __str__(self):
        atts = "\n ".join(["{} ({}) = {} cm".format(attname,attinfo["material"],attinfo["thickness"]) for attname,attinfo in self.attenuators.items()])
        return "Attenuators:\n {}".format(atts)
        
    def addtoconfig(self,config): 
        for attname,attinfo in self.attenuators.items():
            self.addattenuator(config,attname,attinfo["material"],attinfo["thickness"])

    def addattenuator(self,config,attname,material,thickness):
        matname = self.addmaterial(config,material)
        config["attenuators"][attname] = [1,matname,material.density,thickness,1.0]
            
    def addmaterial(self,config,material):
        matname,v = material.pymcaformat()
        if material.nelements>1:
            config["materials"][matname] = v
        return matname
        
    def loadfromconfig(self,config):
        self.attenuators = {}
        for attname,attinfo in self.attenuators["attenuators"].items():
            self.parseattenuator(config,attname,attinfo)

    def parseattenuator(self,config,attname,attinfo):
        matname,density,thickness = attinfo[1],attinfo[2],attinfo[3]
        material = parsematerial(config,matname,density)
        self.attenuators[attname] = {"material":material,"thickness":thickness}
        
    def parsematerial(self,config,matname,density):
        if matname in config["materials"]:
            material = mixture.pymcaformat(config["materials"][matname])
        else:
            material = compoundfromformula.CompoundFromFormula(matname,density)
        return material

    @staticmethod
    def isbefore(name):
        return "BeamFilter" in name

    @staticmethod
    def isafter(name):
        return "BeamFilter" not in name and name != "Detector"
    
    def addtofisx(self,setup,cfg):
        att_after = []
        att_before = []

        for attname,attinfo in self.attenuators.items():
            material = attinfo["material"]
            cfg.addMaterial(material)
            
            if self.isafter(attname):
                att_after.append([material.pymcaname,material.density,attinfo["thickness"],1.0])
            elif self.isbefore(attname):
                att_before.append([material.pymcaname,material.density,attinfo["thickness"],1.0])
        
        if att_before:
            setup.setBeamFilters(att_before)
        if att_after:
            setup.setAttenuators(att_after)
        
    def filter_absorbance(self,energy,source=False):
        if source:
            atts = [self.attenuators[att] for att in self.attenuators if self.isbefore(att)]
        else:
            atts = [self.attenuators[att] for att in self.attenuators if self.isafter(att)]
        if atts:
            return sum([att["material"].mass_att_coeff(energy)*(att["thickness"]*att["material"].density) for att in atts])
        else:
            return np.zeros_like(energy,dtype=float)

    def filter_transmission(self,energy,source=False):
        return np.exp(-self.filter_absorbance(energy,source=source))
        
    def filter_attenuation(self,energy,source=False):
        return 1-self.filter_transmission(energy,source=source)
        
        
class PymcaConfig(object):

    def __init__(self,sample=None,energy=None,flux=None,time=None,fluolinefilter=False):
        self.sample = sample
        self.energy = energy # keV
        self.flux = flux # ph/s
        self.time = time # s

    @property
    def _energy(self):
        return units.magnitude(self.energy,"keV")
    
    @property
    def _time(self):
        return units.magnitude(self.time,"s")
    
    @property
    def _flux(self):
        return units.magnitude(self.flux,"hertz")
        
    def beam(self,config):
        energy = self._energy
        if False:
            config["fit"]["energy"] = [energy, 100*energy]
            config["fit"]["energyflag"] = [1,1-self.fluolinefilter]
            config["fit"]["energyscatter"] = [1,0]
            config["fit"]["energyweight"] = [1e+100,1e-05]
        else:
            config["fit"]["energy"] = [energy,0]
            config["fit"]["energyflag"] = [1,0]
            config["fit"]["energyscatter"] = [1,0]
            config["fit"]["energyweight"] = [1.0,0.]
        config["fit"]["scatterflag"] = 1
        
        config["concentrations"]["flux"] = self._flux
        config["concentrations"]["time"] = self._time
        
    def background(self,config):
        config["fit"]["stripflag"] = 1
        config["fit"]["stripalgorithm"] = 1
        config["fit"]["snipwidth"] = 100
    
    def peakshape(self,config):
        config["fit"]["hypermetflag"] = 3 # shorttail
        config["fit"]["hypermetflag"] = 1

    def fill(self,config):
        self.sample.addtoconfig(config,self._energy)
        self.beam(config)
        self.background(config)
        self.peakshape(config)


class FisxConfig():
    
    FISXMATERIALS = None
    
    def init(self):
        if self.FISXMATERIALS is None:
            self.FISXMATERIALS = fisx.Elements()
            self.FISXMATERIALS.initializeAsPyMca()
        
    @contextlib.contextmanager
    def init_ctx(self):
        self.init()
        yield
        
    def setDetector(self,setup,detector):
        with self.init_ctx():
            detector.addtofisx(setup,self.FISXMATERIALS)

    def addMaterial(self,material):
        if material.nelements>1:
            with self.init_ctx():
                self.FISXMATERIALS.addMaterial(material.tofisx(),errorOnReplace=False)
        return material.pymcaname


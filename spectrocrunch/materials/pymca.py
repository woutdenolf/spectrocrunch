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

    def addMaterial(self,config,material):
        matname,v = material.pymcaformat()
        if material.nelements>1:
            config["materials"][matname] = v
        return matname
        
        
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


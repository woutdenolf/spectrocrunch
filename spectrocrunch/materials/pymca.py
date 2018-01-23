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

    def __init__(self,sample=None,emin=0.9,emax=None,\
                energy=None,weights=None,scatter=None,flux=None,time=None):
        self.sample = sample
        
        self.energy = units.magnitude(energy,"keV")
        self.emin = emin
        self._emax = emax
        self.weights = weights
        self.scatter = scatter
        
        self.flux = units.magnitude(flux,"hertz")
        self.time = units.magnitude(time,"s")
    
    @property
    def emax(self):
        if self._emax is None:
            return np.max(self.energy)
        else:
            return self._emax
        
    def beam(self,config):
        # TODO: move to source?
        config["fit"]["energy"] = np.asarray(self.energy).tolist()
        config["fit"]["energyweight"] = np.asarray(self.weights).tolist()
        config["fit"]["energyscatter"] = np.asarray(self.scatter,dtype=int).tolist()
        config["fit"]["energyflag"] = np.ones_like(self.energy,dtype=int).tolist()
        config["fit"]["scatterflag"] = 1
        
        config["concentrations"]["flux"] = self.flux
        config["concentrations"]["time"] = self.time
        
    def background(self,config):
        config["fit"]["stripflag"] = 1
        config["fit"]["stripalgorithm"] = 1
        config["fit"]["snipwidth"] = 100
    
    def peakshape(self,config):
        config["fit"]["hypermetflag"] = 3 # shorttail
        config["fit"]["hypermetflag"] = 1

    def fill(self,cfg):
        self.sample.addtopymca(self,cfg)
        self.beam(cfg)
        self.background(cfg)
        self.peakshape(cfg)
        print(cfg)

    def addtopymca_material(self,cfg,material):
        matname,v = material.topymca()
        if material.nelements>1:
            cfg["materials"][matname] = v
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

    def addtofisx_material(self,material):
        if material.nelements>1:
            with self.init_ctx():
                self.FISXMATERIALS.addMaterial(material.tofisx(),errorOnReplace=False)
        return material.pymcaname


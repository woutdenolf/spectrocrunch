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

from ..common import instance

import numpy as np
import itertools


class Layer(object):

    def __init__(self,material=None,thickness=None,fixed=False,ml=None):
        """
        Args:
            material(list(spectrocrunch.materials.compound|mixture)): material composition
            thickness(num): thickness in micron
            fixed(bool): thickness and composition are fixed
            ml(Multilayer): part of this ensemble
        """
        self.material = material
        self.thickness = thickness
        self.fixed = fixed
        self.ml = ml

    def __str__(self):
        return "{}: {} μm".format(self.material,self.thickness)
        
    def __getattr__(self,name):
        try:
            return getattr(self.material,name)
        except AttributeError as e:
            return getattr(self.ml,name)
            
    @property
    def xraythickness(self):
        return self.thickness*self.ml.cosanglein*1e-4
    
    @xraythickness.setter
    def xraythickness(self,value):
        self.thickness = value/self.ml.cosanglein*1e-4
    
    def absorbance(self,energy,fine=False,decomposed=False):
        if decomposed:
            return {"cs":self.mass_att_coeff(energy,fine=fine,decomposed=decomposed),"thickness":self.xraythickness,"density":self.density}
        else:
            return self.mass_att_coeff(energy,fine=fine,decomposed=decomposed)*(self.xraythickness*self.density)


class Multilayer(object):
    """
    Class representing a multilayer of compounds or mixtures
    """
    
    def __init__(self, material=None, thickness=None, detector=None):
        """
        Args:
            material(list(spectrocrunch.materials.compound|mixture)): layer composition
            thickness(list(num)): layer thickness in micron
            detector(spectrocrunch.detectors.xrf.Detector): 
        """
        
        self.detector = detector

        if not instance.isarray(material):
            material = [material]
        if not instance.isarray(thickness):
            thickness = [thickness]
        self.layers = [Layer(material=mat,thickness=t,ml=self) for mat,t in itertools.izip(material,thickness)]
    
    def __getattr__(self,name):
        return getattr(self.detector,name)
     
    def __len__(self):
        return len(self.layers)
    
    @property
    def nlayers(self):
        return len()
        
    def __iter__(self):
        return iter(self.layers)
    
    def fixediter(self):
        for layer in self:
            if layer["fixed"]:
                yield layer
    
    def freeiter(self):
        for layer in self:
            if not layer["fixed"]:
                yield layer
            
    def __getitem__(self,index):
        return self.layers[index]
  
    def __str__(self):
        layers = "\n ".join("L{}. {}".format(i,str(layer)) for i,layer in enumerate(self))
        return "Multilayer (ordered top-bottom):\n {}".format(layers)

    def markscatterer(self,name):
        for layer in self:
            layer.markscatterer(name)

    def ummarkscatterer(self):
        for layer in self:
            layer.ummarkscatterer()
   
    def markabsorber(self,symb,shells=[],fluolines=[]):
        """
        Args:
            symb(str): element symbol
        """
        for layer in self:
            layer.markabsorber(symb,shells=shells,fluolines=fluolines)

    def unmarkabsorber(self):
        for layer in self:
            layer.unmarkabsorber()
    
    def absorbance(self,energy,fine=False,decomposed=False):
        if decomposed:
            return [layer.absorbance(energy,fine=fine,decomposed=decomposed) for layer in self]
        else:
            return np.sum([layer.absorbance(energy,fine=fine,decomposed=decomposed) for layer in self],axis=0)
            
    def transmission(self,energy,fine=False,decomposed=False):
        A = self.absorbance(energy,fine=fine,decomposed=decomposed)
        if decomposed:
            return A
        else:
            return np.exp(-A)

    def fixlayers(self,ind=None):
        if ind is None:
            for layer in self:
                layer["fixed"] = True
        else:
            for i in ind:
                self[i]["fixed"] = True

    def freelayers(self,ind=None):
        if ind is None:
            for layer in self:
                layer["fixed"] = False
        else:
            for i in ind:
                self[i]["fixed"] = False

    def refinethickness(self,energy,absorbance,constant=False,constraint=True):
        y = absorbance

        A = [layer.density*layer.mass_att_coeff(energy) for layer in self.freeiter]
        
        for layer in self.fixediter:
            y = y-layer.density*layer.xraythickness*layer.mass_att_coeff(energy)
        
        if constant:
            A.append(np.ones_like(energy))
            A = np.vstack(A).T
            
            if constraint:
                lb = np.zeros(len(A),dtype=float)
                lb[-1] = -np.inf
                ub = np.inf
                thickness = fit1d.lstsq_bound(A,y,lb,ub)
            else:
                thickness = fit1d.lstsq(A,y)
            thickness = thickness[:-1]
        else:
            A = np.vstack(A).T
            if constraint:
                thickness = fit1d.lstsq_nonnegative(A,y)
            else:
                thickness = fit1d.lstsq(A,y)
                
        for t,layer in itertools.izip(thickness,self.freeiter):
            layers.xraythickness = t

    def spectrum(self,energy):
        # TODO: generate true spectrum with the fisx library
        
        # Oversimplified: isotropic scattering of all absorbed radiation
        detfrac = self.detector.solidangle/(4*np.pi)
        genfrac = 1-self.transmission(energy)
        return energy,detfrac*genfrac



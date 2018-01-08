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

import numpy as np
import pandas as pd
import scipy.integrate
import scipy.special

from ..common import instance
from ..common import cache
from ..common import listtools
from . import xrayspectrum

class Layer(object):

    def __init__(self,material=None,thickness=None,fixed=False,ml=None):
        """
        Args:
            material(list(spectrocrunch.materials.compound|mixture)): material composition
            thickness(num): thickness in cm
            fixed(bool): thickness and composition are fixed
            ml(Multilayer): part of this ensemble
        """
        self.material = material
        self.thickness = thickness
        self.fixed = fixed
        self.ml = ml

    def __str__(self):
        return "{} μm ({})".format(self.thickness*1e4,self.material)
        
    def __getattr__(self,attr):
        try:
            return getattr(self.material,attr)
        except AttributeError:
            return getattr(self.ml,attr)
            
    @property
    def xraythickness(self):
        return self.thickness/self.ml.cosnormin
    
    @xraythickness.setter
    def xraythickness(self,value):
        self.thickness = value*self.ml.cosnormin
    
    def absorbance(self,energy,fine=False,decomposed=False):
        if decomposed:
            return {"cs":self.material.mass_att_coeff(energy,fine=fine,decomposed=decomposed),"thickness":self.xraythickness,"density":self.material.density}
        else:
            return self.material.mass_att_coeff(energy,fine=fine,decomposed=decomposed)*(self.xraythickness*self.material.density)

class Multilayer(cache.Cache):
    """
    Class representing a multilayer of compounds or mixtures
    """
    
    def __init__(self, material=None, thickness=None, detector=None):
        """
        Args:
            material(list(spectrocrunch.materials.compound|mixture)): layer composition
            thickness(list(num)): layer thickness in cm
            detector(spectrocrunch.detectors.xrf.Detector): 
        """
        
        self.detector = detector

        if not instance.isarray(material):
            material = [material]
        if not instance.isarray(thickness):
            thickness = [thickness]
        self.layers = [Layer(material=mat,thickness=t,ml=self) for mat,t in zip(material,thickness)]
        
        super(Multilayer,self).__init__(force=True)
    
    def __getattr__(self,name):
        return getattr(self.detector,name)
     
    def __len__(self):
        return len(self.layers)
    
    @property
    def nlayers(self):
        return len(self.layers)
        
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
        layers = "\n ".join("Layer {}. {}".format(i,str(layer)) for i,layer in enumerate(self))
        return "Multilayer (ordered top-bottom):\n {}".format(layers)

    def markscatterer(self,name):
        for layer in self:
            layer.markscatterer(name)

    def ummarkscatterer(self):
        for layer in self:
            layer.ummarkscatterer()
   
    @property
    def density(self):
        return np.vectorize(lambda layer:layer.density)(self)
    
    @property
    def thickness(self):
        return np.vectorize(lambda layer:layer.thickness)(self)
    
    def mass_att_coeff(self,energy):
        return np.asarray([layer.mass_att_coeff(energy) for layer in self])
        
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
                
        for t,layer in zip(thickness,self.freeiter):
            layers.xraythickness = t

    def _cache_layerinfo(self):
        t = np.empty(self.nlayers+1)
        np.cumsum(self.thickness,out=t[1:])
        t[0] = 0
        if self.reflection:
            zexit = 0.
        else:
            zexit = t[-1]
        return {"cumul_thickness":t,"zexit":zexit}

    def _zlayer(self,z):
        """Get layer in which z falls
        
        Args:
            z(num): depth
        
        Returns:
            num: 
                0 when z<=0
                n+1 when z>totalthickness
                {1,...,n} otherwise (the layers)
        """
        layerinfo = self.getcashed("layerinfo")
        return np.asscalar(np.digitize(z, layerinfo["cumul_thickness"], right=True))
    
    def _cache_attenuationinfo(self,energy):
        energy = np.sort(instance.asarray(energy))
        nenergies = len(energy)

        density = self.density[:,np.newaxis]
        thickness = self.thickness[:,np.newaxis]
        mu = self.mass_att_coeff(energy).reshape((self.nlayers,nenergies))
        
        # We will add one layer at the beginning and one at the end, both vacuum
        
        # linear attenuation coefficient for each layer
        linatt = mu * density
        linattout = np.empty((self.nlayers+2,nenergies),dtype=linatt.dtype)
        linattout[1:-1,:] = linatt
        linattout[[0,-1],:] = 0 # outside sample (vacuum)
        
        # Cumulative linear attenuation coefficient (= linatt*z + correction)
        attall = (linatt*thickness).sum(axis=0)
        A = np.empty((self.nlayers+2,nenergies),dtype=attall.dtype)
        A[0,:] = 0 # before sample (vacuum)
        A[-1,:] = attall # after sample
        
        for i in range(nenergies):
            tmp = np.subtract.outer(linatt[:,i],linatt[:,i])
            tmp *= thickness
            A[1:-1,i] = np.triu(tmp).sum(axis=0)

        linattout = pd.DataFrame(linattout,columns=energy,index=range(self.nlayers+2))
        A = pd.DataFrame(A,columns=energy,index=range(self.nlayers+2))
        
        return {"linatt":linattout,"linatt_cumulcor":A}
    
    def _cum_attenuation(self,z,energy):
        """Total attenuation from surface to z
        
        Args:
            z(num): depth of attenuation
            energy(num|array): energies to be attenuation
            
        Returns:
            array:
        """
        
        lz = self._zlayer(z)
        att = self.getcashed("attenuationinfo")

        att = z * att["linatt"].loc[lz][energy] + att["linatt_cumulcor"].loc[lz][energy]
        if np.isscalar(att):
            return np.asscalar(att)
        else:
            return att.values

    def _transmission(self,zi,zj,cosaij,energy):
        """Transmission from depth zi to zj
        
        Args:
            zi(num): start depth of attenuation
            zj(num): end depth of attenuation
            cosaij(num): angle with surface normal
            energy(array): energies to be attenuation
            
        Returns:
            array:
        """

        datt = self._cum_attenuation(zj,energy)-self._cum_attenuation(zi,energy)
        return np.exp(-datt/cosaij)

    def _cache_interactioninfo(self,energy,emin=None,emax=None,scatteringangle=None,ninteractions=None):
    
        probabilities = [None]*ninteractions
        energyindex = [None]*ninteractions

        energy = instance.asarray(energy)
        energy.sort()

        for i in range(ninteractions):
            def f(x,energy=energy):
                return (np.abs(energy-x)).argmin()
            
            interactions = [layer.xrayspectrum(energy,emin=emin,emax=emax) for layer in self]
            probs = [dict(interaction.probabilities) for interaction in interactions]
            
            interactions = list(set(listtools.flatten(p.keys() for p in probs)))
            energy = list(set(listtools.flatten(interaction.energy(scatteringangle=scatteringangle) for interaction in interactions)))
            energy.sort()

            probabilities[i] = probs
            energyindex[i] = f

        return {"probabilities":probabilities,"energyindex":energyindex,"energies":energy}

    def _gentransmission(self,zi,zj,cosaij,i,energyi,energyj,interactionj):
        """Generation at depth zi and then transmission from zi to zj
        
        Args:
            zi(num): start depth of attenuation
            zj(num): end depth of attenuation
            cosaij(num): angle with surface normal
            i(num): interaction 1, 2, ...
            line(): energies to be attenuation
            
        Returns:
            array:
        """
        lz = self._zlayer(zi)
        if lz==0:
            return 0
        
        interactions = self.getcashed("interactioninfo")

        try:
            prob = interactions["probabilities"][i-1][lz-1][interactionj]
        except:
            return 0
        
        ind = interactions["energyindex"][i-1](energyi)
        return prob[ind]*self._transmission(zi,zj,cosaij,energyj)

    def _gentransmission_saintegrated(self,zi,zj,i,energyi,energyj,interactionj):
        """Generation at depth zi and then transmission from zi to zj
        
        Args:
            zi(num): start depth of attenuation
            zj(num): end depth of attenuation
            aij(num): angle with surface normal
            i(num): interaction 1, 2, ...
            line(): energies to be attenuation
            
        Returns:
            array:
        """
        lz = self._zlayer(zi)
        if lz==0:
            return 0
        
        interactions = self.getcashed("interactioninfo")

        try:
            prob = interactions["probabilities"][i-1][lz-1][interactionj]
        except:
            return 0
        
        ind = interactions["energyindex"][i-1](energyi)

        datt = self._cum_attenuation(zj,energyj)-self._cum_attenuation(zi,energyj)
        return prob[ind]*scipy.special.exp1(datt)*(2*np.pi)

    @cache.withcache("layerinfo")
    def xrayspectrum(self,energy,emin=0,emax=None):

        spectrum = xrayspectrum.Spectrum()

        a01 = self.detector.anglenormin
        a12 = self.detector.anglenormout
        scatteringangle1 = a12-a01
        cosafirst = np.cos(a01)
        cosalast = np.cos(a12)
        layerinfo = self.getcashed("layerinfo")
        za = layerinfo["cumul_thickness"][0]
        zb = layerinfo["cumul_thickness"][-1]
        zfirst = layerinfo["cumul_thickness"][0]
        zlast = layerinfo["zexit"]
        
        def addsources(data):
            data = instance.asarray(data)
            while len(data.shape)>=2:
                data = data.sum(axis=-1)
            return data*(self.detector.solidangle/cosafirst)

        with self.cachectx("interactioninfo",energy,emin=emin,emax=emax,ninteractions=2,scatteringangle=scatteringangle1):
            interactions = self.getcashed("interactioninfo")
            energies = interactions["energies"]

            getinteractions = lambda x: list(set(listtools.flatten(p.keys() for p in x)))


            with self.cachectx("attenuationinfo",energies):

                interactions1 = getinteractions(interactions["probabilities"][0])
                gen1 = {}
                path = lambda z1: self._transmission(zfirst,z1,cosafirst,energy0)*\
                                  self._gentransmission(z1,zlast,cosalast,1,energy0,energy1,interaction1)

                for energy0 in instance.asarray(energy):
                    for interaction1 in interactions1:
                        energy1 = interaction1.energy(scatteringangle=scatteringangle1)
                        
                        gen = [scipy.integrate.quad(path, za, zb)[0] for energy1 in instance.asarray(energy1)]

                        gen1[interaction1] = addsources(gen)
                print(gen1)

                if True:
                    interactions2 = getinteractions(interactions["probabilities"][1])
                    gen2 = {}
                    path = lambda z1,z2: self._transmission(zfirst,z1,cosafirst,energy0)*\
                                         self._gentransmission_saintegrated(z1,z2,1,energy0,energy1,interaction1)*\
                                         self._gentransmission(z2,zlast,cosalast,2,energy1,energy2,interaction2)

                    for energy0 in instance.asarray(energy):
                        for interaction2 in interactions2:
                            scatteringangle2 = scatteringangle1
                            energy2 = interaction2.energy(scatteringangle=scatteringangle2)
                            for interaction1 in interactions1:
                                energy1 = interaction1.energy(scatteringangle=scatteringangle1)

                                gen = [[scipy.integrate.nquad(path, [(za, zb)]*2)[0] for energy1 in instance.asarray(energy1)] for energy2 in instance.asarray(energy2)]
                                
                                if interaction2 in gen2:
                                    gen2[interaction2] += addsources(gen)
                                else:
                                    gen2[interaction2] = addsources(gen)

            
                    print(gen2)

        spectrum.cs = gen1
        spectrum.xlim = [min(energies),max(energies)]
        spectrum.density = 1
        spectrum.title = str(self)
        spectrum.type = spectrum.TYPES.interactionyield

        return spectrum





            

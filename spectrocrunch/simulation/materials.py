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

from .simul import with_simulmetaclass

import numpy as np
import contextlib

from . import noisepropagation
from ..materials.compoundfromlist import compoundfromlist as compound
from ..materials.mixture import mixture
from ..materials.types import fractionType
from ..math.fit1d import lstsq
from ..common.Enum import Enum
from ..common.instance import isarray

from ..common import listtools

interactionType = Enum(["transmission","fluorescence","elastic","inelastix"])

class Multilayer(with_simulmetaclass()):
    """
    Class representing an area material
    """
    
    def __init__(self, material=None, thickness=None, anglein=None, angleout=None):
        """
        Args:
            material(list(spectrocrunch.materials.compound|mixture)): material composition
            thickness(num): thickness in micron
            anglein(num): angle (deg) between primary beam and surface normal (pointing inwards)
            angleout(num): angle (deg) between fluorescene path to detector and surface normal (pointing inwards)
        """
        self.required(material,"material")
        self.required(thickness,"thickness")
        self.required(anglein,"anglein")
        self.required(angleout,"angleout")
        
        if hasattr(thickness,"__iter__"):
            self.thickness = np.asarray(thickness,dtype=float)
        else:
            self.thickness = np.asarray([thickness],dtype=float)
        self.thickness *= 1e-4
        
        if hasattr(material,"__iter__"):
            self.material = material
        else:
            self.material = [material]
            
        self.cosanglein = np.cos(np.radians(anglein))
        self.cosangleout = np.cos(np.radians(angleout))
        
        self.nlayers = len(self.material)
        
        self._cache = {}
        
    def layerproperties(self,prop):
        return np.asarray([getattr(self.material[layer],prop) for layer in range(self.nlayers)])
            
    def layercs(self,method,energy,**kwargs):
        cs = [getattr(self.material[layer],method)(energy,**kwargs) for layer in range(self.nlayers)]
        if "decomposed" not in kwargs:
            cs = np.asarray(cs)
        return cs
        
    def layermethod(self,method,*args,**kwargs):
        for layer in range(self.nlayers):
            getattr(self.material[layer],method)(*args,**kwargs)
    
    def markscatterer(self,name):
        self.layermethod("markscatterer",name)
    
    def ummarkscatterer(self):
        self.layermethod("ummarkscatterer")
   
    def markabsorber(self,symb,shells=[],fluolines=[]):
        """
        Args:
            symb(str): element symbol
        """
        self.layermethod("markabsorber",symb,shells=shells,fluolines=fluolines)

    def unmarkabsorber(self):
        self.layermethod("unmarkabsorber")

    def prob_elastic(self,energy):
        density = self.layerproperties("density")
        thickness = self.thickness*self.cosanglein
        mu = self.layercs("mass_att_coeff",energy)
        muscat = self.layercs("rayleigh_cross_section",energy)

        # PATCH: select specific scatterer
        cs = self.material[1].rayleigh_cross_section(energy,decomposed=True)
        muscat = 0.
        for k in cs:
            if k=="verdigris":
                w = cs[k]["w"]
                for e,v in cs[k]["elements"].items():
                    muscat += w*v["w"]*v["cs"]
        muscat = np.asarray([0,muscat,0])
        
        if mu.ndim==2:
            thickness = thickness.reshape((self.nlayers,1))
            density = density.reshape((self.nlayers,1))

        att = mu*density*thickness
        Tin = np.exp(-att/self.cosanglein)
        Tout = np.exp(-att/self.cosangleout)
        
        #chi = mu*(1./self.cosanglein-1./self.cosangleout)
        #player = muscat*(1-np.exp(-chi*density*thickness))/self.cosangleout
        #ind = chi!=0
        #player[ind] /= chi[ind]
        #print player
        #if self.cosangleout>0:
        #    player *= np.exp(-att/self.cosangleout)
        
        player = muscat*density*thickness*np.exp(-att)
        
        print Tin
        print Tout
        for i in range(self.nlayers):
            player[i] *= np.prod(Tin[:i-1,...],axis=0)*np.prod(Tout[i+1:,...],axis=0)
        
        return sum(player)
        
    def prob_transmission(self,energy):
        """Transmission probability for one layer
        """
        T = [np.exp(-self.material[layer].density*self.thickness[layer]*self.cosanglein*self.material[layer].mass_att_coeff(energy)) for layer in range(self.nlayers)]

        #density = self.layerproperties("density")
        #thickness = self.thickness*self.cosanglein
        #mu = self.layercs("mass_att_coeff",energy)
        #if mu.ndim==2:
        #   thickness = thickness.reshape((self.nlayers,1))
        #   density = density.reshape((self.nlayers,1))
        #T = np.exp(-mu*density*thickness)
        
        return np.prod(T,axis=0)
        
    def refinethickness(self,energy,absorbance,layerfixed=None):
        if layerfixed is None:
            layerfixed = []
        y = absorbance

        A = [self.material[layer].density*self.material[layer].mass_att_coeff(energy) for layer in range(self.nlayers) if layer not in layerfixed]

        for layer in range(self.nlayers):
            if layer in layerfixed:
                y -= self.material[layer].density*self.material[layer].mass_att_coeff(energy)

        A = np.vstack(A).T
        thickness = lstsq(A,y)
        ind = [layer not in layerfixed for layer in range(self.nlayers)]
        self.thickness[ind] = thickness

    def absorbance(self,energy):
        return sum([self.material[layer].density*self.thickness[layer]*self.material[layer].mass_att_coeff(energy) for layer in range(self.nlayers)])

    def propagate(self,N,energy,interaction=interactionType.transmission,withnoise=True,forward=True):
        """Error propagation of a number of photons.
               
        Args:
            N(num|array): incomming number of photons with uncertainties
            energy(num|array): energies

        Returns:
            num|numpy.array
        """
        # Bernouilli processes: compounding is the same as multiplication
        #                       so we can multiply the probabilities
        if interaction==interactionType.transmission:
            probsuccess = self.prob_transmission(energy)
        elif interaction==interactionType.fluorescence:
            raise RuntimeError("Fluorescence not implemented yet")
        elif interaction==interactionType.elastic:
            probsuccess = self.prob_elastic(energy)
        elif interaction==interactionType.inelastic:
            raise RuntimeError("Inelastic scattering not implemented yet")

        N,probsuccess = self.propagate_broadcast(N,probsuccess)

        if withnoise:
            process = noisepropagation.bernouilli(probsuccess)
            Nout = noisepropagation.compound(N,process,forward=forward)
        else:
            if forward:
                Nout = N*probsuccess
            else:
                Nout = N/probsuccess
                
        return Nout

    def __str__(self):
        layers = "\n ".join("{}. {}: {} μm".format(i,m,t) for i,(m,t) in enumerate(zip(self.material,self.thickness)))
        return "Multilayer (ordered top-bottom):\n {}".format(layers)

    @contextlib.contextmanager
    def _layer_cache_ctx(self,energy):
        # Already cashed or not
        cached = False
        if "layers" in self._cache:
            cached = self._cache["layers"]["cached"]
        else:
            self._cache["layers"] = {"cached":cached}
        
        # Cache
        if not cached:
            self._layers_cache(energy)
            self._cache["layers"]["cached"] = True
        
        # Use cache
        yield
        
        # Keep cache?
        self._cache["layers"]["cached"] = cached
        
    def _layers_cache(self,energy):
        """
        
        Args:
            energy(num|array): energies to be attenuation
            
        Returns:
            None
        """

        density = self.layerproperties("density")
        thickness = self.thickness
        mu = self.layercs("mass_att_coeff",energy)
        
        # nlayers x nenergy
        if mu.ndim!=2:
            mu = mu.reshape((self.nlayers,1))
        nlayers,nenergies = mu.shape
        thickness = thickness.reshape((nlayers,1))
        density = density.reshape((nlayers,1))
        
        # total attenuation
        murho = mu * density
        attall = (murho*thickness).sum(axis=0)
        att = np.empty((nlayers+2,nenergies),dtype=attall.dtype)
        att[0,:] = 0
        att[-1,:] = attall
        
        # for each layer: the attenuation of previous layers 
        #                 minus the attenuation of previous layers
        #                 if they all had the same composition as this layer
        for i in range(nenergies):
            tmp = np.subtract.outer(murho[:,i],murho[:,i])
            tmp *= thickness
            att[1:-1,i] = np.triu(tmp).sum(axis=0)
        
        # nlayers+2 x nenergy
        tmp = np.empty((nlayers+2,nenergies),dtype=attall.dtype)
        tmp[1:-1,:] = murho
        tmp[[0,-1],:] = 0
        self._cache["layers"]["murho"] = tmp   # as if all previous layers had the same composition
        self._cache["layers"]["cor"] = att     # correction on this assumption
    
        # borders
        t = np.empty(self.nlayers+1)
        np.cumsum(self.thickness,out=t[1:])
        t[0] = 0
        self._cache["layers"]["t"] = t
    
    def _zlayer(self,z):
        """
        Args:
            z(num|array): depth
        
        Returns:
            num|array: 
                0 when z<=0
                n+1 when z>totalthickness
                {1,...,n} otherwise
        """
        return np.digitize(z, self._cache["layers"]["t"], right=True)
        
    def _cum_attenuation(self,z,energy):
        """
        Args:
            z(num|array): depth of attenuation
            energy(num|array): energies to be attenuation
            
        Returns:
            array: nz x nenergy
        """
        
        with self._layer_cache_ctx(energy):
            layers = self._zlayer(z)
            att = z.reshape((z.size,1)) * self._cache["layers"]["murho"][layers,:] + self._cache["layers"]["cor"][layers,:]
            return att
        
    def _transmission(self,zi,zj,cosaij,energy):
        """
        Args:
            zi(num|array): start depth of attenuation (dims: nz)
            zj(num|array): end depth of attenuation (dims: nz)
            cosaij(num|array): angle with surface normal (dims: nz)
            energy(num|array): energies to be attenuation (dims: nenergy)
            
        Returns:
            array: nz x nenergy
        """
        
        with self._layer_cache_ctx(energy):
            datt = self._cum_attenuation(zj,energy)-self._cum_attenuation(zi,energy)
            if isarray(cosaij) and datt.ndim==2:
                cosaij = cosaij.reshape((cosaij.size,1))
            return np.exp(-datt/cosaij)
        
    def _fluorescence(self,energy):
        mufluo = self.layercs("fluorescence_cross_section",energy,decomposed=True)
        
        interactions = [{}]*self.nlayers
        for l in range(self.nlayers):
            for k,v in mufluo[l].items():
                if "elements" in v:
                    for k,v in v["elements"].items():
                        if k not in lines:
                            if k.isabsorber():
                                lines[k] = list(listtools.flatten(s.fluolines for s in k.shells))
                else:
                    if k not in lines:
                        if k.isabsorber():
                                lines[k] = list(listtools.flatten(s.fluolines for s in k.shells))
                    
        print lines    
    
classes = Multilayer.clsregistry
aliases = Multilayer.aliasregistry
factory = Multilayer.factory



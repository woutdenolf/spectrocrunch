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
import itertools

from . import noisepropagation
from ..materials.compoundfromlist import CompoundFromList as compound
from ..materials.mixture import Mixture
from ..materials.multilayer import Multilayer as _Multilayer
from ..materials.types import fractionType
from ..materials import element
from ..materials import interaction
from ..math import fit1d
from ..common.Enum import Enum
from ..common.instance import isarray
from ..common import listtools

interactionType = Enum(["transmission","fluorescence","elastic","inelastix"])

class Multilayer(with_simulmetaclass(_Multilayer)):
    """
    Class representing an area material
    """
    
    def __init__(self, **kwargs):
        super(Multilayer, self).__init__(**kwargs)
        
        self._cache = {}

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
        
    def propagate(self,N,energy,interaction=interactionType.transmission,forward=True):
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
            probsuccess = self.transmission(energy)
        elif interaction==interactionType.fluorescence:
            raise RuntimeError("Fluorescence not implemented yet")
        elif interaction==interactionType.elastic:
            probsuccess = self.prob_elastic(energy)
        elif interaction==interactionType.inelastic:
            raise RuntimeError("Inelastic scattering not implemented yet")

        N,probsuccess = self.propagate_broadcast(N,probsuccess)

        if noisepropagation.israndomvariable(N):
            process = noisepropagation.bernouilli(probsuccess)
            Nout = noisepropagation.compound(N,process,forward=forward)
        else:
            if forward:
                Nout = N*probsuccess
            else:
                Nout = N/probsuccess
                
        return Nout

    @contextlib.contextmanager
    def _interaction_cache_ctx(self,topic,*args,**kwargs):
        # Already cashed or not
        cached = False
        if "interaction" in self._cache:
            cached = self._cache["interaction"]["cached"]
        else:
            self._cache["interaction"] = {"cached":cached}
        
        # Cache
        if not cached:
            if topic=="attenuation":
                self._attenuation_cache(*args,**kwargs)
            elif topic=="generation":
                self._generation_cache(*args,**kwargs)
                
            self._cache["interaction"]["cached"] = True
        
        # Use cache
        yield
        
        # Keep cache?
        self._cache["interaction"]["cached"] = cached

    def _generation_cache(self,energy):
        """
        Args:
            energy(num|array): primary beam energy
            
        Returns:
            None
        """
        
        if not isarray(energy):
            energy = [energy]
        
        finteractions = set()
        for mat in self.material:
            finteractions = finteractions | set(listtools.flatten(mat.fluointeractions()))   
        
        source = set(interaction.InteractionSource(e,i) for i,e in enumerate(energy))
        
        
        n = 3
        
        for i in range(1,n+1):
            interactions = finteractions
            interactions = interactions | set(interaction.InteractionElScat(s) for s in source)
            interactions = interactions | set(interaction.InteractionInelScat(s,45) for s in source) # TODO: angle not fixed!
            
            source = sorted(interactions)
            print source
            print np.array([l.energy for l in source])

        #self._attenuation_cache(energy)
        
        #print self.layercs("fluorescence_cross_section",energy,decomposed=True)
        
        #print self.layercs("rayleigh_cross_section",energy,decomposed=True)
        
        #print self.layercs("compton_cross_section",energy,decomposed=True)
        
    def _attenuation_cache(self,energy):
        """
        Args:
            energy(num|array): energies to be attenuation
            
        Returns:
            None
        """

        # Borders
        t = np.empty(self.nlayers+1)
        np.cumsum(self.thickness,out=t[1:])
        t[0] = 0
        self._cache["interaction"]["t"] = t
        
        # Attenuation
        linatt,cor = self._cumulated_linear_attenuation_coefficient(energy)
        self._cache["interaction"]["linatt"] = linatt # nlayers+2 x nenergy
        self._cache["interaction"]["linatt_cumulcor"] = cor # nlayers+2 x nenergy
        
    def _cumulated_linear_attenuation_coefficient(self,energy):
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
        
        # linear attenuation coefficient for each layer
        linatt = mu * density
        linattout = np.empty((nlayers+2,nenergies),dtype=linatt.dtype)
        linattout[1:-1,:] = linatt
        linattout[[0,-1],:] = 0 # outside sample (assume no attenuation)

        # cumulative linear attenuation coefficient (= linatt*z + correction)
        attall = (linatt*thickness).sum(axis=0)
        A = np.empty((nlayers+2,nenergies),dtype=attall.dtype)
        A[0,:] = 0
        A[-1,:] = attall
        
        for i in range(nenergies):
            tmp = np.subtract.outer(linatt[:,i],linatt[:,i])
            tmp *= thickness
            A[1:-1,i] = np.triu(tmp).sum(axis=0)

        return linattout,A
    
    def _zlayer(self,z):
        """Get layer in which z falls
        
        Args:
            z(num|array): depth
        
        Returns:
            num|array: 
                0 when z<=0
                n+1 when z>totalthickness
                {1,...,n} otherwise
        """
        return np.digitize(z, self._cache["interaction"]["t"], right=True)
        
    def _cum_attenuation(self,z,energy):
        """Total attenuation from surface to z
        
        Args:
            z(num|array): depth of attenuation
            energy(num|array): energies to be attenuation
            
        Returns:
            array: nz x nenergy
        """
        
        with self._interaction_cache_ctx("attenuation",energy):
            lz = self._zlayer(z)
            att = z.reshape((z.size,1)) * self._cache["interaction"]["linatt"][lz,:] + self._cache["interaction"]["linatt_cumulcor"][lz,:]
            return att
        
    def _transmission(self,zi,zj,cosaij,energy):
        """Transmission from depth zi to zj
        
        Args:
            zi(num|array): start depth of attenuation (dims: nz)
            zj(num|array): end depth of attenuation (dims: nz)
            cosaij(num|array): angle with surface normal (dims: nz)
            energy(num|array): energies to be attenuation (dims: nenergy)
            
        Returns:
            array: nz x nenergy
        """
        
        with self._interaction_cache_ctx("attenuation",energy):
            datt = self._cum_attenuation(zj,energy)-self._cum_attenuation(zi,energy)
            if isarray(cosaij) and datt.ndim==2:
                cosaij = cosaij.reshape((cosaij.size,1))
            return np.exp(-datt/cosaij)

    def _fluorescence(self,energy):
        """
        
        Args:
            energy(num|array): primary beam energy
        """
        
        
        with self._interaction_cache_ctx("generation",energy):
            pass
                
        
        

    
classes = Multilayer.clsregistry
aliases = Multilayer.aliasregistry
factory = Multilayer.factory



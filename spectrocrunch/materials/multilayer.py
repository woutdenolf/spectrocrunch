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
import collections
import fisx

from ..common import instance
from ..common import cache
from ..common import listtools
from . import xrayspectrum
from ..simulation.classfactory import with_metaclass
from ..simulation import noisepropagation
from . import pymca
from . import element

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

    def __getattr__(self,attr):
        return getattr(self.material,attr)

    def __str__(self):
        return "{} μm ({})".format(self.thickness*1e4,self.material)
  
    @property
    def xraythickness(self):
        return self.thickness/self.ml.geometry.cosnormin
    
    @xraythickness.setter
    def xraythickness(self,value):
        self.thickness = value*self.ml.geometry.cosnormin
    
    def absorbance(self,energy,fine=False,decomposed=False):
        if decomposed:
            return {"cs":self.material.mass_att_coeff(energy,fine=fine,decomposed=decomposed),\
                    "thickness":self.xraythickness,\
                    "density":self.density}
        else:
            return self.material.mass_att_coeff(energy,fine=fine,decomposed=decomposed)*\
                    (self.xraythickness*self.density)

    def absorbanceout(self,energy):
        return self.material.mass_att_coeff(energy)*\
                    (self.thickness/self.ml.cosnormout*self.density)

    def addtofisx(self,setup,cfg):
        name = cfg.addtofisx_material(self.material)
        return [name,self.density,self.thickness]

    def fisxgroups(self,emin=0,emax=np.inf):
        return self.material.fisxgroups(emin=emin,emax=emax)
    
    
class Multilayer(with_metaclass(cache.Cache)):
    """
    Class representing a multilayer of compounds or mixtures
    """
    
    FISXCFG = pymca.FisxConfig()

    def __init__(self, material=None, thickness=None, geometry=None):
        """
        Args:
            material(list(spectrocrunch.materials.compound|mixture)): layer composition
            thickness(list(num)): layer thickness in cm
            geometry(spectrocrunch.geometries.base.PointGeometry): 
        """
        
        self.geometry = geometry

        if not instance.isarray(material):
            material = [material]
        if not instance.isarray(thickness):
            thickness = [thickness]
        self.layers = [Layer(material=mat,thickness=t,ml=self) for mat,t in zip(material,thickness)]
        
        super(Multilayer,self).__init__(force=True)
    
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
    
    @property
    def xraythickness(self):
        return np.vectorize(lambda layer:layer.xraythickness)(self)
        
    def mass_att_coeff(self,energy):
        """Total mass attenuation coefficient
        
        Args:
            energy(num|array): keV
        Returns:
            array: nz x nenergy
        """
        return np.asarray([instance.asarray(layer.mass_att_coeff(energy)) for layer in self])
        
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
   
    def absorbanceout(self,energy):
        return np.sum([layer.absorbanceout(energy) for layer in self],axis=0)
                     
    def transmission(self,energy,fine=False,decomposed=False):
        A = self.absorbance(energy,fine=fine,decomposed=decomposed)
        if decomposed:
            return A
        else:
            return np.exp(-A)

    def transmissionout(self,energy):
        return np.exp(-self.absorbanceout(energy))
        
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
        if self.geometry.reflection:
            zexit = 0.
        else:
            zexit = t[-1]
        return {"cumul_thickness":t,"zexit":zexit}

    def _zlayer(self,z):
        """Get layer in which z falls
        
        Args:
            z(num|array): depth
        
        Returns:
            num|array: 
                0 when z<=0
                n+1 when z>totalthickness
                {1,...,n} otherwise (the layers)
        """
        layerinfo = self.getcashed("layerinfo")
        ret = np.digitize(z, layerinfo["cumul_thickness"], right=True)
        try:
            return np.asscalar(ret)
        except:
            return ret
    
    def _cache_attenuationinfo(self,energy):
        energy = np.unique(instance.asarray(energy))
        nenergies = len(energy)

        density = self.density[:,np.newaxis]
        thickness = self.thickness[:,np.newaxis]
        mu = self.mass_att_coeff(energy)
        
        # We will add one layer at the beginning and one at the end, both vacuum
        
        # linear attenuation coefficient for each layer
        linatt = mu * density
        linattout = np.empty((self.nlayers+2,nenergies),dtype=linatt.dtype)
        linattout[1:-1,:] = linatt
        linattout[[0,-1],:] = 0 # outside sample (vacuum)
        
        # Cumulative linear attenuation coefficient (= linatt*z + correction)
        attall = (linatt*thickness).sum(axis=0)
        cor = np.empty((self.nlayers+2,nenergies),dtype=attall.dtype)
        cor[0,:] = 0 # before sample (vacuum)
        cor[-1,:] = attall # after sample
        
        for i in range(nenergies):
            tmp = np.subtract.outer(linatt[:,i],linatt[:,i])
            tmp *= thickness
            cor[1:-1,i] = np.triu(tmp).sum(axis=0)

        linattout = pd.DataFrame(linattout,columns=energy,index=range(self.nlayers+2))
        cor = pd.DataFrame(cor,columns=energy,index=range(self.nlayers+2))
        
        return {"linatt":linattout,"linatt_cumulcor":cor}
    
    def _cum_attenuation(self,z,energy):
        """Total attenuation from surface to z
        
        Args:
            z(num|array): depth of attenuation
            energy(num|array): energies to be attenuation
            
        Returns:
            array: nz x nenergy
        """
        
        lz = self._zlayer(z)
        att = self.getcashed("attenuationinfo")
        linatt = att["linatt"].loc[lz][energy]
        cor = att["linatt_cumulcor"].loc[lz][energy]
        if linatt.ndim!=0:
            linatt = linatt.values
            cor = cor.values
        if linatt.ndim==2:
            z = z[:,np.newaxis]
        
        return z*linatt + cor

    def _transmission(self,zi,zj,cosaij,energy):
        """Transmission from depth zi to zj
        
        Args:
            zi(num|array): start depth of attenuation (nz)
            zj(num|array): end depth of attenuation (nz)
            cosaij(num|array): angle with surface normal (nz)
            energy(num|array): energies to be attenuation (nenergy)
            
        Returns:
            array: nz x nenergy
        """
        datt = self._cum_attenuation(zj,energy)-self._cum_attenuation(zi,energy)
        if datt.ndim==2:
            if instance.isarray(cosaij):
                cosaij = cosaij[:,np.newaxis]
        #assert(sum(instance.asarray(-datt/cosaij)>0)==0)
        return np.exp(-datt/cosaij)

    def _cache_interactioninfo(self,energy,emin=None,emax=None,ninteractions=None,geomkwargs=None):
        # Pepare resulting lists
        _nlayers = self.nlayers+2
        _ninteractions = ninteractions+2
        probabilities = [None]*_ninteractions
        energyindex = [None]*_ninteractions

        # Source energies
        energy,func = instance.asarrayf(energy)
        nsource = len(energy)
        probabilities[0] = pd.DataFrame(columns=[xrayspectrum.RayleighLine(energy)])
        
        # Calculate interaction probabilities (ph/cm/srad)
        def getenergy(x,**kwargs):
            return list(listtools.flatten(interaction.energy(**kwargs) for interaction in x.columns))
        
        for i in range(ninteractions):
            energy = getenergy(probabilities[i],**geomkwargs)
            nenergy = len(energy)
            
            def f(x,energy=energy):
                return (np.abs(energy-x)).argmin()

            # Interaction probabilities:
            #  column -> interaction
            #  index -> [layer,source]
            probs = [None]*_nlayers
            probs[1:-1] = [pd.DataFrame.from_items(layer.xrayspectrum(energy,emin=emin,emax=emax).probabilities) for layer in self]
            probs[0] = pd.DataFrame(index=range(nenergy))
            probs[-1] = probs[0]
            probs = pd.concat(probs)
            probs.fillna(0., inplace=True)
            probs.index = pd.MultiIndex.from_product([np.arange(_nlayers), range(nenergy)],names=["layer","source"])

            probabilities[i+1] = probs
            energyindex[i+1] = f

        return {"probabilities":probabilities,\
                "energyindex":energyindex,\
                "getenergy":getenergy}
        
    def _genprobabilities(self,zi,i,energyi,interactionj):
        """Generation at depth zi
        
        Args:
            zi(num|array): start depth of attenuation
            i(num): interaction 1, 2, ...
            energyi(num): energy in
            interactionj(object|array): 
            
        Returns:
            array:
        """
        
        lz = self._zlayer(zi)
        lzarr = instance.isarray(lz)
        if lzarr:
            lz = lz.tolist()
            
        interactioninfo = self.getcashed("interactioninfo")
        energyindex = interactioninfo["energyindex"][i](energyi)
        
        # Advanced indexing on MultiIndex: does not preserve order and repeats
        probs = interactioninfo["probabilities"][i].loc[(lz,energyindex),interactionj]
        if probs.ndim!=0:
            if lzarr:
                # apply order and repeats in lz (assume energyindex is a scalar)
                probs.index = probs.index.droplevel(1)
                probs = probs.loc[lz]
            probs = probs.values
            
        return probs
            
    def _gentransmission(self,zi,zj,cosaij,i,energyi,energyj,interactionj):
        """Generation at depth zi and then transmission from zi to zj
        
        Args:
            zi(num|array): start depth of attenuation
            zj(num|array): end depth of attenuation
            cosaij(num|array): angle with surface normal
            i(num): interaction 1, 2, ...
            energyi(num): energy in
            energyj(num|array): energy out
            interactionj(object|array): 
            
        Returns:
            array:
        """
        probs = self._genprobabilities(zi,i,energyi,interactionj)
        T = self._transmission(zi,zj,cosaij,energyj)
        return probs*T

    def _gentransmission_saintegrated(self,zi,zj,i,energyi,energyj,interactionj):
        """Generation at depth zi, then transmission from zi to zj and then integrate
           over the solid angle of emission from zi to zj (hemisphere).
        
        Args:
            zi(num|array): start depth of attenuation
            zj(num|array): end depth of attenuation
            i(num): interaction 1, 2, ...
            energyi(num): energy in
            energyj(num): energy out
            interactionj(): energies to be attenuation
            
        Returns:
            array:
        """
        
        # assume isotropic radiation
        probs = self._genprobabilities(zi,i,energyi,interactionj)
        Aj = self._cum_attenuation(zj,energyj)
        Ai = self._cum_attenuation(zi,energyj)

        barri = instance.isarray(zi)
        barrj = instance.isarray(zj)
        if barri and barrj:
            probs = instance.asarray(probs)[:,np.newaxis]
            Ai = instance.asarray(Ai)[:,np.newaxis]
            Aj = instance.asarray(Aj)[np.newaxis,:]
            
        #func = lambda theta,phi: probs*np.exp(-(Aj-Ai)/np.cos(theta))*np.tan(theta)
        #return np.nquad(func,[(0,np.pi/2),(0,2*np.pi)])
        
        return (2*np.pi)*probs*scipy.special.exp1(Aj-Ai)
                
    def _primary_interaction(self,selfabs=True):
        """Returns the ph/srad generated per source line after 1 interaction (without efficiency term)
        """
        interactionindex = 1

        interactioninfo = self.getcashed("interactioninfo")
        energy0 = interactioninfo["getenergy"](interactioninfo["probabilities"][0])
        
        nsource = len(energy0)
        interactions1 = interactioninfo["probabilities"][interactionindex].columns
        nlayers = self.nlayers
        nlines = len(interactions1)
        
        # Integrated attenuation over the sample thickness
        if selfabs:
            geomkwargs = self.geometry.xrayspectrumkwargs()
            energy1 = interactioninfo["getenergy"](interactioninfo["probabilities"][interactionindex],**geomkwargs)
            
            att = self.getcashed("attenuationinfo")
            cosafirst = self.geometry.cosnormin
            cosalast = self.geometry.cosnormout
            mu0 = att["linatt"][energy0].values/cosafirst
            mu1 = att["linatt"][energy1].values/cosalast
            cor0 = att["linatt_cumulcor"][energy0].values/cosafirst
            cor1 = att["linatt_cumulcor"][energy1].values/cosalast

            chi = mu1[1:-1,:,np.newaxis] - mu0[1:-1,np.newaxis,:]
            chicor = cor1[1:-1,:,np.newaxis] - cor0[1:-1,np.newaxis,:]
            
            layerinfo = self.getcashed("layerinfo")
            J2 = np.exp(chi*layerinfo["cumul_thickness"][1:,np.newaxis,np.newaxis])
            J2 -= np.exp(chi*layerinfo["cumul_thickness"][:-1,np.newaxis,np.newaxis])
            J2 /= chi
            J2 *= np.exp(chicor)
            if not self.geometry.reflection:
                J2 *= np.exp(-cor1[-1,np.newaxis,:,np.newaxis])
            
            # nlayers x nenergy1 x nenergy0 -> nlayers x nsource x nenergy1
            J2 = np.transpose(J2,[0,2,1])
            
            # nlayers x nenergy0 x nenergy1 -> nlayers x nsource x nlines (reduce scattering lines)
            interactions1exp = list(listtools.flatten([interaction]*interaction.nenergy for interaction in interactions1))
            indC = np.asarray([interaction=="Compton" for interaction in interactions1exp])
            indR = np.asarray([interaction=="Rayleigh" for interaction in interactions1exp])
            indF = ~indC & ~indR
            
            indsource = range(nsource)
            J2 = np.concatenate((J2[...,indF],\
                               J2[:,indsource,indC][...,np.newaxis],\
                               J2[:,indsource,indR][...,np.newaxis]),axis=-1)
            interactions1 = interactions1.tolist()
            interactions1.insert(nlines, interactions1.pop(interactions1.index("Compton")))
            interactions1.insert(nlines, interactions1.pop(interactions1.index("Rayleigh")))
        else:
            # nlayers x 1 x 1
            J2 = self.thickness[:,np.newaxis,np.newaxis]

        # Interaction probability: nlayers x nsource x nlines  (ph/cm/srad)
        probs = interactioninfo["probabilities"][interactionindex].loc[(range(1,self.nlayers+1),),interactions1]
        probs = probs.values.reshape((nlayers,nsource,nlines))

        # Rate: fluoresence/scattering per incoming photon
        J2 = J2*probs # (ph/cm/srad) -> (ph/srad)
        
        # Sum over layers
        J2 = J2.sum(axis=0).T # nlines x nsource
        
        # Sum over source lines for fluorescence
        J2 = dict(zip(interactions1,J2))

        return J2
        
    def _primary_interaction_numerical(self):
        """Returns the ph/srad generated per source line after 1 interaction (without efficiency term)
        """
        interactionindex = 1
        
        cosafirst = self.geometry.cosnormin
        cosalast = self.geometry.cosnormout
        layerinfo = self.getcashed("layerinfo")
        za = layerinfo["cumul_thickness"][0]
        zb = layerinfo["cumul_thickness"][-1]
        zfirst = layerinfo["cumul_thickness"][0]
        zlast = layerinfo["zexit"]
        geomkwargs = self.geometry.xrayspectrumkwargs()
        
        interactioninfo = self.getcashed("interactioninfo")
        energy0 = interactioninfo["getenergy"](interactioninfo["probabilities"][0])
        interactions1 = interactioninfo["probabilities"][interactionindex].columns
        
        J2 = {}
        path = lambda z1: self._transmission(zfirst,z1,cosafirst,en0)*\
                       self._gentransmission(z1,zlast,cosalast,interactionindex,en0,en1,interaction1)
        
        def numintegrate(path,za,zb):
            return scipy.integrate.quad(path, za, zb)[0]
        
        n = (zb-za)/min(self.thickness)*100
        def numintegratefast(path,za,zb):
            x = np.linspace(za,zb,n)
            y = path(x)
            return np.trapz(y, x=x)
            #return scipy.integrate.trapz(y, x=x)
            
        #import matplotlib.pyplot as plt

        for interaction1 in interactions1:
            energy1 = interaction1.energy(**geomkwargs)
            if isinstance(interaction1,xrayspectrum.FluoZLine):
                energy1 = [energy1]*len(energy0)

            gen = [numintegrate(path, za, zb)\
                    for en0,en1 in zip(energy0,energy1)]
            
            #plt.figure()
            #x = np.linspace(za,zb,n)
            #plt.plot(x,path(x))
            #plt.show()
            
            J2[interaction1] = np.asarray(gen)

        return J2
    
    def _secondary_interaction_numerical(self):
        """Returns the ph/srad generated per source line after 2 interactions (without efficiency term)
        """
        interactionindex = 2
        
        cosafirst = self.geometry.cosnormin
        cosalast = self.geometry.cosnormout
        layerinfo = self.getcashed("layerinfo")
        za = layerinfo["cumul_thickness"][0]
        zb = layerinfo["cumul_thickness"][-1]
        zfirst = layerinfo["cumul_thickness"][0]
        zlast = layerinfo["zexit"]
        geomkwargs1 = self.geometry.xrayspectrumkwargs()
        geomkwargs2 = geomkwargs1
        
        interactioninfo = self.getcashed("interactioninfo")
        energy0 = interactioninfo["getenergy"](interactioninfo["probabilities"][0])
        interactions1 = interactioninfo["probabilities"][1].columns
        interactions2 = interactioninfo["probabilities"][2].columns
        J3 = {}
        path = lambda z1,z2: self._transmission(zfirst,z1,cosafirst,en0)[:,np.newaxis]*\
                             self._gentransmission_saintegrated(z1,z2,interactionindex-1,en0,en1,interaction1)*\
                             self._gentransmission(z2,zlast,cosalast,interactionindex,en1,en2,interaction2)[np.newaxis,:]

        def numintegrate(path,za,zb):
            return scipy.integrate.nquad(path, [(za, zb)]*2)[0]
        
        n = (zb-za)/min(self.thickness)*100
        def numintegratefast(path,za,zb):
            x1 = np.linspace(za,zb,n)
            x2 = np.linspace(za,zb,n)
            y = path(x1,x2)
            y = np.trapz(y, x=x1, axis=0)
            y = np.trapz(y, x=x2, axis=0)
            return y
        
        import matplotlib.pyplot as plt
        
        for interaction1 in interactions1:
            energy1 = interaction1.energy(**geomkwargs1)
            if isinstance(interaction1,xrayspectrum.FluoZLine):
                energy1 = [energy1]*len(energy0)

            for interaction2 in interactions2:
                energy2 = interaction2.energy(**geomkwargs2)
                if isinstance(interaction2,xrayspectrum.FluoZLine):
                    energy2 = [energy2]*len(energy1)

                for en0,en1,en2 in zip(energy0,energy1,energy2):
                    x1 = np.linspace(za,zb,n)
                    x2 = np.linspace(za,zb,n+1)
                    
                    print(self._transmission(zfirst,x1,cosafirst,en0).shape)
                    print(self._gentransmission_saintegrated(x1,x2,interactionindex-1,en0,en1,interaction1).shape)
                    print(self._gentransmission(x2,zlast,cosalast,interactionindex,en1,en2,interaction2).shape)
                    
                    plt.figure()
                    
                    img = path(x1,x2)
                    plt.imshow(img)
                    plt.show()
                    
                gen = [numintegrate(path, za, zb)\
                    for en0,en1,en2 in zip(energy0,energy1,energy2)]

                J3[interaction2] = np.asarray(gen)
  
        return J3

    def addtofisx(self,setup,cfg):
        cfg.init()
        setup.setSample([layer.addtofisx(setup,cfg) for layer in self])
        self.geometry.addtofisx(setup,cfg)

    def addtopymca_matrix(self,setup,cfg,name):
        anglein = self.geometry.anglein
        angleout = self.geometry.angleout 
        scatteringangle = self.geometry.scatteringangle 
        if name=="MULTILAYER":
            density,thickness = 0.,0.
        else:
            v = cfg["materials"][name]
            density,thickness = v["Density"], v["Thickness"]
        cfg["attenuators"]["Matrix"] = [1, name, density, thickness, anglein, angleout, 0, scatteringangle]
    
    def addtopymca_layer(self,setup,cfg,index,layer):
        name = setup.addtopymca_material(cfg,layer)
        l = "Layer{}".format(index)
        cfg["multilayer"][l] = [1,name,layer.density,layer.thickness]

    def addtopymca_shells(self,setup,cfg,elements):
        emax = setup.emax
        emin = setup.emin
        
        for e in elements:
            shells = xrayspectrum.Shell.pymcafactory((e.Z,emin,emax))
            if shells:
                cfg["peaks"][str(e)] = shells

    def addtopymca(self,setup,cfg):
        if self.nlayers==1:
            name = setup.addtopymca_material(cfg,self[0])
            self.addtopymca_shells(setup,cfg,self[0].elements)
            self.addtopymca_matrix(setup,cfg,name)
        else:
            for index,layer in enumerate(self):
                self.addlayer(setup,cfg,index,layer)
                self.addtopymca_shells(setup,cfg,layer.elements)
            self.addtopymca_matrix(setup,cfg,'MULTILAYER')
        self.geometry.addtopymca(setup,cfg)
    
    def _parse_fisx_result(self,gen):
        result = {}
        
        # Get fluorescence rates from fisx (add escape peaks)
        for group in gen:
            el = element.Element(group.split(' ')[0])
            for layer in gen[group]:
                for peak in gen[group][layer]:
                    line = xrayspectrum.FluoLine(peak.split(' ')[0])
                    line = xrayspectrum.FluoZLine(el,line)
                    rate = gen[group][layer][peak]["rate"]
                    if line in result:
                        result[line] += rate
                    else:
                        result[line] = rate
        
        # Correction for detector in transmission
        if not self.geometry.reflection:
            for line in result:
                energy = line.energy(**self.geometry.xrayspectrumkwargs())
                result[line] = result[line]*self.transmissionout(energy)
        
        return result

    def _dict_to_spectrum(self,gen,emin=0,emax=None):
        if emax is None:
            allenergies = list(listtools.flatten(line.energy(**self.geometry.xrayspectrumkwargs()) for line in gen))  
            emax = max(allenergies)
            
        spectrum = xrayspectrum.Spectrum()
        spectrum.update(gen)
        spectrum.xlim = [emin,emax]
        spectrum.density = None
        spectrum.title = str(self)
        spectrum.type = spectrum.TYPES.interactionyield
        return spectrum
        
    def _print_fisx(self,fluo):
        print("Element   Peak          Energy       Rate      Secondary      Tertiary      Efficiency")
        for key in fluo:
            for layer in fluo[key]:
                peakList = list(fluo[key][layer].keys())
                peakList.sort()
                for peak in peakList:
                    # energy of the peak
                    energy = fluo[key][layer][peak]["energy"]
                    # expected measured rate
                    rate = fluo[key][layer][peak]["rate"]
                    # primary photons (no attenuation and no detector considered)
                    primary = fluo[key][layer][peak]["primary"]
                    # secondary photons (no attenuation and no detector considered)
                    secondary = fluo[key][layer][peak]["secondary"]
                    # tertiary photons (no attenuation and no detector considered)
                    tertiary = fluo[key][layer][peak].get("tertiary", 0.0)
                    efficiency = fluo[key][layer][peak].get("efficiency", 0.0)
                    # correction due to secondary excitation
                    enhancement2 = (primary + secondary) / primary
                    enhancement3 = (primary + secondary + tertiary) / primary
                    print("%s   %s    %.4f     %.3g     %.5g    %.5g    %.5g" % \
                                       (key, peak + (13 - len(peak)) * " ", energy,
                                       rate, enhancement2, enhancement3, efficiency))
    
    @staticmethod
    def _parse_weights(weights,nsource):
        if weights is None:
            if nsource==1:
                weights = 1.
            else:
                weights = np.ones(nsource,dtype=float)/nsource
        else:
            weights,func = instance.asarrayf(weights,dtype=float)
            weights = func(weights/weights.sum())

        return weights

    def _interactions_fisx(self,energy0,weights,ninteractions,emin=0,emax=None):
        energy0 = instance.asarray(energy0)
        nsource = len(energy0)
        if emax is None:
            emax = np.max(energy0)
            
        setup = fisx.XRF()
        
        # Add sample, detector and geometry
        self.addtofisx(setup,self.FISXCFG)
        
        # Add source
        setup.setBeam(energy0,weights=self._parse_weights(weights,nsource))

        # Peak groups
        groups = {}
        for layer in self:
            groups.update(layer.fisxgroups(emin=emin,emax=emax)) 
        
        def shellparse(shell):
            shell = str(shell)
            if shell.startswith("M"):
                shell = "M"
            return shell
            
        groups = set(["{} {}".format(el,shellparse(shell)) for el,shells in groups.items() for shell in shells])

        # Get fluorescence
        secondary = 2*(ninteractions>1)
        gen = setup.getMultilayerFluorescence(\
                        groups,self.FISXCFG.FISXMATERIALS,\
                        secondary=secondary,useMassFractions=1)
                        
        #self._print_fisx(gen)
        result = self._parse_fisx_result(gen)
        
        return result

    def _interactions_applydetector(self,gen):
        """Convert ph/srad to ph and apply detector efficiency
        """
        
        geom = self.geometry.xrayspectrumkwargs()
        
        lines = gen.keys()
        energysource = lines[lines.index("Rayleigh")].energy(**geom)
        energydet = [k.energy(**geom) for k in lines]
        ind = np.cumsum([listtools.length(en) for en in energydet])
        ind = np.insert(ind,0,0)
        ind = zip(ind[:-1],ind[1:])

        energydet = list(listtools.flatten(energydet))
        efficiency = self.geometry.efficiency(energysource,energydet)

        for k,(a,b) in zip(lines,ind):
            if a+1==b: # Fluorescence
                eff = efficiency[:,a]
            else: # Scattering
                eff = np.diag(efficiency[:,a:b])
            gen[k] = gen[k]*eff
            
    @cache.withcache("layerinfo")
    def xrayspectrum(self, energy0, emin=0, emax=None, method="analytical", ninteractions=1, weights=None, scattering=True):
    
        if method=="fisx":
            if scattering:
                method="analytical"
            #if not self.geometry.reflection and self.nlayers>1:
            #    method="analytical"
            if not self.geometry.reflection:
                method="analytical"
            
            if ninteractions>3:
                method="analytical"
                
        if method=="analytical":
            if ninteractions>=2:
                method="numerical"
    
        if method=="fisx":
            gen = self._interactions_fisx(energy0,weights,ninteractions,emin=emin,emax=emax)
        else:
            geomkwargs=self.geometry.xrayspectrumkwargs()
            with self.cachectx("interactioninfo",energy0,emin=emin,emax=emax,\
                                ninteractions=ninteractions,\
                                geomkwargs=geomkwargs):
                interactioninfo = self.getcashed("interactioninfo")
                allenergies = interactioninfo["getenergy"](interactioninfo["probabilities"][-2],**geomkwargs)
                with self.cachectx("attenuationinfo",allenergies):
                    # Primary interaction (with self-absorption)
                    if method=="numerical":
                        gen = self._primary_interaction_numerical()
                    else:
                        gen = self._primary_interaction()
                    
                    # Secondary interaction (with self-absorption)
                    if ninteractions>=2:
                        for k,v in self._secondary_interaction_numerical().items():
                            if k in gen:
                                gen[k] += v
                            else:
                                gen[k] = v

            # Apply filter attenuation (source and detection) + detector efficiency
            self._interactions_applydetector(gen)

        spectrum = self._dict_to_spectrum(gen,emin=emin,emax=emax)

        if method!="fisx":
            spectrum.apply_weights(weights)
            #spectrum.sum_sources() # this is already done when fisx is used
            
        return spectrum

    def propagate(self,N,energy,interaction="transmission",forward=True):
        """Error propagation of a number of photons.
               
        Args:
            N(num|array): incomming number of photons with uncertainties
            energy(num|array): energies

        Returns:
            num|numpy.array
        """
        # Bernouilli processes: compounding is the same as multiplication
        #                       so we can multiply the probabilities
        if interaction=="transmission":
            probsuccess = self.transmission(energy)
        else:
            raise RuntimeError("{} not implemented yet".format(interaction))

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



        
                 
classes = Multilayer.clsregistry
aliases = Multilayer.aliasregistry
factory = Multilayer.factory


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
            return {"cs":self.material.mass_att_coeff(energy,fine=fine,decomposed=decomposed),\
                    "thickness":self.xraythickness,\
                    "density":self.material.density}
        else:
            return self.material.mass_att_coeff(energy,fine=fine,decomposed=decomposed)*\
                    (self.xraythickness*self.material.density)

    def addtofisx(self,setup,cfg):
        name = cfg.addMaterial(self.material)
        return [name,self.material.density,self.thickness]

    def fisxgroups(self,emin=0,emax=np.inf):
        return self.material.fisxgroups(emin=emin,emax=emax)
    
    
class Multilayer(with_metaclass(cache.Cache)):
    """
    Class representing a multilayer of compounds or mixtures
    """
    
    FISXCFG = pymca.FisxConfig()

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
        return np.exp(-datt/cosaij)

    def _cache_interactioninfo(self,energy,emin=None,emax=None,scatteringangle=None,ninteractions=None):
        # Pepare resulting lists
        _nlayers = self.nlayers+2
        _ninteractions = ninteractions+2
        probabilities = [None]*_ninteractions
        energyindex = [None]*_ninteractions

        # Source energies
        energy = instance.asarray(energy)
        probabilities[0] = pd.DataFrame(columns=[xrayspectrum.RayleighLine(energy)])
        
        # Calculate interaction probabilities (ph/cm/srad)
        getenergy = lambda x: list(listtools.flatten(interaction.energy(scatteringangle=scatteringangle) for interaction in x.columns))
        for i in range(ninteractions):
            energy = getenergy(probabilities[i])
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

        return {"probabilities":probabilities,"energyindex":energyindex,"getenergy":getenergy}

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
            
        return probs*self._transmission(zi,zj,cosaij,energyj)

    def _gentransmission_saintegrated(self,zi,zj,i,energyi,energyj,interactionj):
        """Generation at depth zi and then transmission from zi to zj
        
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
        lz = self._zlayer(zi)
        
        interactioninfo = self.getcashed("interactioninfo")

        ind = interactioninfo["energyindex"][i](energyi)
        
        if instance.isarray(lz):
            probs = np.asarray([interactioninfo["probabilities"][i][k][interactionj][ind] for k in lz])
        else:
            probs = interactioninfo["probabilities"][i][lz][interactionj][ind]

        datt = self._cum_attenuation(zj,energyj)-self._cum_attenuation(zi,energyj)
        
        #func = lambda theta,phi: probs[ind]*np.exp(-datt/np.cos(theta))*np.tan(theta)
        #return np.nquad(func,[(0,np.pi/2),(0,2*np.pi)])
        
        return probs*scipy.special.exp1(datt)*(2*np.pi)

    def _primary_interaction(self,selfabs=True):
        """Returns the ph/srad generated per source line after 1 interaction (without efficiency term)
        """
        interactionindex = 1

        interactioninfo = self.getcashed("interactioninfo")
        energy0 = interactioninfo["getenergy"](interactioninfo["probabilities"][interactionindex-1])
        nsource = len(energy0)
        interactions1 = interactioninfo["probabilities"][interactionindex].columns
        nlayers = self.nlayers
        nlines = len(interactions1)
        
        # Integrated attenuation over the sample thickness
        if selfabs:
            energy1 = interactioninfo["getenergy"](interactioninfo["probabilities"][interactionindex])
            
            att = self.getcashed("attenuationinfo")
            cosafirst = self.detector.cosnormin
            cosalast = self.detector.cosnormout
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
        probs = probs.swaplevel(0,1)
        probs.sort_index(inplace=True)
        probs = probs.values.reshape((nlayers,nsource,nlines))
        
        # Rate: fluoresence/scattering per incoming photon
        J2 = J2*probs # (ph/cm/srad) -> (ph/srad)
        return dict(zip(interactions1,J2.sum(axis=0).T))
        
    def _primary_interaction_numerical(self):
        """Returns the ph/srad generated per source line after 1 interaction (without efficiency term)
        """
        interactionindex = 1
        
        cosafirst = self.detector.cosnormin
        cosalast = self.detector.cosnormout
        layerinfo = self.getcashed("layerinfo")
        za = layerinfo["cumul_thickness"][0]
        zb = layerinfo["cumul_thickness"][-1]
        zfirst = layerinfo["cumul_thickness"][0]
        zlast = layerinfo["zexit"]
        scatteringangle = self.detector.scatteringangle
        
        interactioninfo = self.getcashed("interactioninfo")
        energy0 = interactioninfo["getenergy"](interactioninfo["probabilities"][interactionindex-1])
        interactions1 = interactioninfo["probabilities"][interactionindex].columns
        
        J2 = {}
        path = lambda z1: self._transmission(zfirst,z1,cosafirst,en0)*\
                       self._gentransmission(z1,zlast,cosalast,interactionindex,en0,en1,interaction1)
        
        def numintegrate(path,za,zb):
            return scipy.integrate.quad(path, za, zb)[0]
        
        def numintegratefast(path,za,zb):
            x = np.linspace(za,zb,n)
            y = path(x)
            return scipy.integrate.trapz(y, x=x)
            
        #import matplotlib.pyplot as plt
        #n = (zb-za)/min(self.thickness)*100
        
        for interaction1 in interactions1:
            energy1 = interaction1.energy(scatteringangle=scatteringangle)
            
            if isinstance(interaction1,xrayspectrum.FluoZLine):
                en1 = energy1
                gen = [numintegrate(path, za, zb)\
                        for en0 in energy0]
            else:
                energy1 = instance.asarray(energy1)
                gen = [numintegrate(path, za, zb)\
                        for en0,en1 in zip(energy0,energy1)]
            
            #plt.figure()
            #x = np.linspace(za,zb,n)
            #y = path(x)
            #print(interaction1,gen,numintegratefast(path, za, zb))
            #plt.plot(x,y)
            #plt.show()
            
            J2[interaction1] = np.asarray(gen)
            
        return J2
    
    def _secondary_interaction_numerical(self):
        """Returns the ph/srad generated per source line after 2 interactions (without efficiency term)
        """
        #TODO: not finished
        interactionindex = 2
        
        cosafirst = self.detector.cosnormin
        cosalast = self.detector.cosnormout
        layerinfo = self.getcashed("layerinfo")
        za = layerinfo["cumul_thickness"][0]
        zb = layerinfo["cumul_thickness"][-1]
        zfirst = layerinfo["cumul_thickness"][0]
        zlast = layerinfo["zexit"]
        
        interactioninfo = self.getcashed("interactioninfo")
        interactions2  = interactioninfo["probabilities"][interactionindex].columns
        J3 = {}
        path = lambda z1,z2: self._transmission(zfirst,z1,cosafirst,en0)*\
                             self._gentransmission_saintegrated(z1,z2,interactionindex-1,en0,en1,interaction1)*\
                             self._gentransmission(z2,zlast,cosalast,interactionindex,en1,n2,interaction2)

        def addsources(data):
            data = instance.asarray(data)
            while data.ndim>=2:
                data = data.sum(axis=-1)
            return data

        for energy0 in instance.asarray(energy):
            for interaction2 in interactions2:
                scatteringangle2 = scatteringangle1
                energy2 = interaction2.energy(scatteringangle=scatteringangle2)
                for interaction1 in interactions1:
                    energy1 = interaction1.energy(scatteringangle=scatteringangle1)

                    #plt.imshow(arr)
                    #plt.show()
                    
                    gen = [[[scipy.integrate.nquad(path, [(za, zb)]*2)[0]\
                            for en0 in instance.asarray(energy0)]\
                            for en1 in instance.asarray(energy1)]\
                            for en2 in instance.asarray(energy2)]
                    
                    if interaction2 in gen2:
                        J3[interaction2] += addsources(gen)
                    else:
                        J3[interaction2] = addsources(gen)
                        
        return J3

    def addtofisx(self,setup,cfg):
        cfg.init()
        setup.setSample([layer.addtofisx(setup,cfg) for layer in self])
        self.detector.addtofisx(setup,cfg)

    @staticmethod
    def _parse_fisx_result(gen,ret,nsource,sourceindex):
        for group in gen:
            el = element.Element(group.split(' ')[0])
            for layer in gen[group]:
                for peak in gen[group][layer]:
                    line = xrayspectrum.FluoLine(peak.split(' ')[0])
                    line = xrayspectrum.FluoZLine(el,line)
                    rate = gen[group][layer][peak]["rate"]
                    energy = gen[group][layer][peak]["rate"]
                    if line not in ret:
                        ret[line] = np.zeros(nsource,dtype=type(rate))
                    ret[line][sourceindex] += rate

        return ret
        
    def _dict_to_spectrum(self,gen):
        allenergies = list(listtools.flatten(line.energy(scatteringangle=self.detector.scatteringangle) for line in gen))    
        spectrum = xrayspectrum.Spectrum()
        spectrum.cs = gen
        spectrum.xlim = [min(allenergies),max(allenergies)]
        spectrum.density = 1
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
        
    def _interactions_fisx(self,energy0,ninteractions,emin=0,emax=None):
        energy0 = instance.asarray(energy0)
        nsource = len(energy0)
        if emax is None:
            emax = np.max(energy0)
            
        setup = fisx.XRF()
        self.addtofisx(setup,self.FISXCFG)
        #secondary = 2*(ninteractions>1)
        secondary = 0
        
        groups = {}
        for layer in self:
            groups.update(layer.fisxgroups(emin=emin,emax=emax)) 
        groups = ["{} {}".format(el,shell) for el,shells in groups.items() for shell in shells]
        
        ret = {}
        
        for sourceindex,en0 in enumerate(energy0):
            setup.setBeam(en0)
            gen = setup.getMultilayerFluorescence(\
                        groups,self.FISXCFG.FISXMATERIALS,\
                         secondary=secondary,useMassFractions=1)
            #self._print_fisx(gen)
            self._parse_fisx_result(gen,ret,nsource,sourceindex)

        return ret

    def _interactions_applydetector(self,gen):
        # ph/srad -> ph
        
        lines = gen.keys()
        energysource = lines[lines.index("Rayleigh")].energy
        energy = np.asarray([k.energy(scatteringangle=self.detector.scatteringangle) for k in gen])
        efficiency = self.detector.efficiency(energy,energysource)

        for k,eff in zip(gen,efficiency):
            gen[k] = gen[k]*eff

    @cache.withcache("layerinfo")
    def xrayspectrum(self, energy0, emin=0, emax=None, calc="fisx", ninteractions=1):
    
        if calc=="fisx" and ninteractions<=3:
            gen = self._interactions_fisx(energy0,ninteractions,emin=emin,emax=emax) 
        else:
            with self.cachectx("interactioninfo",energy0,emin=emin,emax=emax,\
                                ninteractions=ninteractions,\
                                scatteringangle=self.detector.scatteringangle):
                interactioninfo = self.getcashed("interactioninfo")
                allenergies = interactioninfo["getenergy"](interactioninfo["probabilities"][-2])
                
                with self.cachectx("attenuationinfo",allenergies):
                    if calc=="numerical":
                        gen = self._primary_interaction_numerical()
                    else:
                        gen = self._primary_interaction()

            self._interactions_applydetector(gen)

        spectrum = self._dict_to_spectrum(gen)

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


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
import logging

from ..utils import instance
from ..utils import cache
from ..utils import listtools
from ..math import fit1d
from ..math.utils import weightedsum
from . import xrayspectrum
from ..simulation.classfactory import with_metaclass
from ..math import noisepropagation
from . import pymca
from . import element
from ..materials import compoundfromdb

logger = logging.getLogger(__name__)


class Layer(object):

    def __init__(self, material=None, thickness=None, fixed=False, ml=None):
        """
        Args:
            material(compound|mixture|str): material composition
            thickness(num): thickness in cm
            fixed(bool): thickness and composition are fixed
            ml(Multilayer): part of this ensemble
        """
        if instance.isstring(material):
            ret = compoundfromdb.factory(material)
            if ret is None:
                raise RuntimeError("Invalid material {}".format(material))
            material = ret
        self.material = material
        self.thickness = thickness
        self.fixed = fixed
        self.ml = ml

    def __getstate__(self):
        return {'material': self.material,
                'thickness': self.thickness,
                'fixed': self.fixed}

    def __setstate__(self, state):
        self.material = state['material']
        self.thickness = state['thickness']
        self.fixed = state['fixed']

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.material == other.material and \
                self.thickness == other.thickness and \
                self.fixed == other.fixed
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __getattr__(self, attr):
        return getattr(self.material, attr)

    def __str__(self):
        return "{} um ({})".format(self.thickness*1e4, self.material)

    @property
    def xraythicknessin(self):
        return self.thickness/self.ml.geometry.cosnormin

    @xraythicknessin.setter
    def xraythicknessin(self, value):
        self.thickness = value*self.ml.geometry.cosnormin

    @property
    def xraythicknessout(self):
        return self.thickness/self.ml.geometry.cosnormout

    @xraythicknessout.setter
    def xraythicknessout(self, value):
        self.thickness = value*self.ml.geometry.cosnormout

    def absorbance(self, energy, weights=None, out=False, **kwargs):
        kwargs.pop("decomposed", None)
        if out:
            thickness = self.xraythicknessout
        else:
            thickness = self.xraythicknessin
        muL = self.material.mass_att_coeff(
            energy, **kwargs)*(thickness*self.density)
        if weights is None:
            return muL
        else:
            return weightedsum(muL, weights=weights)

    def addtofisx(self, setup, cfg):
        name = cfg.addtofisx_material(self.material)
        return [name, self.density, self.thickness]

    def fisxgroups(self, emin=0, emax=np.inf):
        return self.material.fisxgroups(emin=emin, emax=emax)

    def arealdensity(self):
        wfrac = self.material.elemental_massfractions()
        m = self.density*self.thickness
        return dict(zip(wfrac.keys(),
                        np.asarray(list(wfrac.values()))*m))


class Multilayer(with_metaclass(cache.Cache)):
    """
    Class representing a multilayer of compounds or mixtures
    """

    FISXCFG = pymca.FisxConfig()

    def __init__(self, material=None, thickness=None, fixed=False, geometry=None):
        """
        Args:
            material(list(spectrocrunch.materials.compound|mixture)): layer composition
            thickness(list(num)): layer thickness in cm
            fixed(list(num)): do not change this layer
            geometry(spectrocrunch.geometries.base.Centric): 
        """

        self.geometry = geometry

        if not instance.isarray(material):
            material = [material]
        if not instance.isarray(thickness):
            thickness = [thickness]
        if not instance.isarray(fixed):
            fixed = [fixed]
        if len(fixed) != len(material) and len(fixed) == 1:
            fixed = fixed*len(material)

        self.layers = [Layer(material=mat, thickness=t, fixed=f, ml=self)
                       for mat, t, f in zip(material, thickness, fixed)]

        super(Multilayer, self).__init__(force=True)

    def __getstate__(self):
        return {'layers': self.layers, 'geometry': self.geometry}

    def __setstate__(self, state):
        self.layers = state['layers']
        for layer in self.layers:
            layer.ml = self
        self.geometry = state['geometry']

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.layers == other.layers and \
                self.geometry == other.geometry
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        return len(self.layers)

    def __getitem__(self, index):
        return self.layers[index]

    @property
    def nlayers(self):
        return len(self.layers)

    def fixediter(self):
        for layer in self:
            if layer.fixed:
                yield layer

    def freeiter(self):
        for layer in self:
            if not layer.fixed:
                yield layer

    def __str__(self):
        layers = "\n ".join("Layer {}. {}".format(i, str(layer))
                            for i, layer in enumerate(self))
        return "Multilayer (ordered top-bottom):\n {}".format(layers)

    def markscatterer(self, name):
        for layer in self:
            layer.markscatterer(name)

    def ummarkscatterer(self):
        for layer in self:
            layer.ummarkscatterer()

    @property
    def density(self):
        return np.vectorize(lambda layer: layer.density)(self)

    @property
    def thickness(self):
        return np.vectorize(lambda layer: layer.thickness)(self)

    @property
    def xraythicknessin(self):
        return np.vectorize(lambda layer: layer.xraythicknessin)(self)

    @property
    def xraythicknessout(self):
        return np.vectorize(lambda layer: layer.xraythicknessin)(self)

    def arealdensity(self):
        ret = collections.Counter()
        for layer in self:
            ret.update(layer.arealdensity())
        return dict(ret)

    def elemental_massfractions(self):
        ret = self.arealdensity()
        s = sum(ret.values())
        return {el: w/s for el, w in ret.items()}

    def mass_att_coeff(self, energy):
        """Total mass attenuation coefficient

        Args:
            energy(num|array): keV
        Returns:
            array: nz x nenergy
        """
        return np.asarray([instance.asarray(layer.mass_att_coeff(energy)) for layer in self])

    def markabsorber(self, symb, shells=[], fluolines=[]):
        """
        Args:
            symb(str): element symbol
        """
        for layer in self:
            layer.markabsorber(symb, shells=shells, fluolines=fluolines)

    def unmarkabsorber(self):
        for layer in self:
            layer.unmarkabsorber()

    def absorbance(self, energy, weights=None, out=False, fine=False, decomposed=False):
        if decomposed:
            return [layer.absorbance(energy, weights=weights, out=out, fine=fine) for layer in self]
        else:
            return np.sum([layer.absorbance(energy, weights=weights, out=out, fine=fine) for layer in self], axis=0)

    def transmission(self, energy, weights=None, out=False, fine=False, decomposed=False):
        A = self.absorbance(energy, weights=weights, out=out,
                            fine=fine, decomposed=decomposed)
        if decomposed:
            return A
        else:
            return np.exp(-A)

    def fixlayers(self, ind=None):
        if ind is None:
            for layer in self:
                layer.fixed = True
        else:
            for i in ind:
                self[i].fixed = True

    def freelayers(self, ind=None):
        if ind is None:
            for layer in self:
                layer.fixed = False
        else:
            for i in ind:
                self[i].fixed = False

    def _refine_linear(self, A, y, constant=False, constraint=True):
        y = instance.asarray(y)

        if y.size == 1 and len(A) == 1:
            return y/A[0]

        if constant:
            A.append(np.ones_like(y))
            A = np.vstack(A).T

            if constraint:
                lb = np.zeros(len(A), dtype=float)
                lb[-1] = -np.inf
                ub = np.inf
                params = fit1d.lstsq_bound(A, y, lb, ub)
            else:
                params = fit1d.lstsq(A, y)
            params = params[:-1]
        else:
            A = np.vstack(A).T
            if constraint:
                params = fit1d.lstsq_nonnegative(A, y)
            else:
                params = fit1d.lstsq(A, y)

        return params

    def _refinerhod(self, energy, absorbance, refinedattr, fixedattr, weights=None, **kwargs):
        y = absorbance
        for layer in self.fixediter():
            y = y-layer.absorbance(energy)

        A = [layer.mass_att_coeff(energy) for layer in self.freeiter()]
        if weights is not None:
            A = [weightedsum(csi, weights=weights) for csi in A]

        params = self._refine_linear(A, y, **kwargs)
        for param, layer in zip(params, self.freeiter()):
            setattr(layer, refinedattr, param/getattr(layer, fixedattr))
            logger.info('Refined {} of "{}": {}'.format(
                refinedattr, layer, getattr(layer, refinedattr)))

    def refinecomposition(self, energy, absorbance, weights=None, fixthickness=True, **kwargs):
        y = absorbance
        for layer in self.fixediter():
            y = y-layer.absorbance(energy, weights=weights)

        A = []
        for layer in self.freeiter():
            mu = layer.mass_att_coeff(energy, decomposed=True)
            w, cs = layer.csdict_parse(mu)
            if weights is not None:
                cs = [weightedsum(csi, weights=weights) for csi in cs]
            A.extend(cs)

        params = self._refine_linear(A, y, **kwargs)

        for layer in self.freeiter():
            n = layer.nparts
            w = params[0:n]
            params = params[n:]

            s = w.sum()
            w = w/s
            w = dict(zip(layer.parts.keys(), w))
            layer.change_fractions(w, "mass")

            if fixthickness:
                layer.density = s/layer.xraythicknessin
                logger.info(
                    'Refined density of "{}": {} g/cm^3'.format(layer, layer.density))
            else:
                layer.xraythicknessin = s/layer.density
                logger.info(
                    'Refined thickness "{}": {} g/cm^3'.format(layer, layer.xraythicknessin))

    def refinethickness(self, energy, absorbance, **kwargs):
        self._refinerhod(energy, absorbance,
                         "xraythicknessin", "density", **kwargs)

    def refinedensity(self, energy, absorbance, **kwargs):
        self._refinerhod(energy, absorbance, "density",
                         "xraythicknessin", **kwargs)

    def _cache_layerinfo(self):
        t = np.empty(self.nlayers+1)
        np.cumsum(self.thickness, out=t[1:])
        t[0] = 0
        if self.geometry.reflection:
            zexit = 0.
        else:
            zexit = t[-1]
        return {"cumul_thickness": t, "zexit": zexit}

    def _zlayer(self, z):
        """Get layer in which z falls

        Args:
            z(num|array): depth

        Returns:
            num|array: 
                0 when z<=0
                n+1 when z>totalthickness
                {1,...,n} otherwise (the layers)
        """
        layerinfo = self.getcache("layerinfo")
        ret = np.digitize(z, layerinfo["cumul_thickness"], right=True)
        return instance.asscalar(ret)

    def _cache_attenuationinfo(self, energy):
        energy = np.unique(instance.asarray(energy))
        nenergies = len(energy)

        density = self.density[:, np.newaxis]
        thickness = self.thickness[:, np.newaxis]
        mu = self.mass_att_coeff(energy)

        # We will add one layer at the beginning and one at the end, both vacuum

        # linear attenuation coefficient for each layer
        linatt = mu * density
        linattout = np.empty((self.nlayers+2, nenergies), dtype=linatt.dtype)
        linattout[1:-1, :] = linatt
        linattout[[0, -1], :] = 0  # outside sample (vacuum)

        # Cumulative linear attenuation coefficient (= linatt*z + correction)
        attall = (linatt*thickness).sum(axis=0)
        cor = np.empty((self.nlayers+2, nenergies), dtype=attall.dtype)
        cor[0, :] = 0  # before sample (vacuum)
        cor[-1, :] = attall  # after sample

        for i in range(nenergies):
            tmp = np.subtract.outer(linatt[:, i], linatt[:, i])
            tmp *= thickness
            cor[1:-1, i] = np.triu(tmp).sum(axis=0)

        linattout = pd.DataFrame(
            linattout, columns=energy, index=range(self.nlayers+2))
        cor = pd.DataFrame(cor, columns=energy, index=range(self.nlayers+2))

        return {"linatt": linattout, "linatt_cumulcor": cor}

    def _cum_attenuation(self, z, energy):
        """Total attenuation from surface to z

        Args:
            z(num|array): depth of attenuation
            energy(num|array): energies to be attenuation

        Returns:
            array: nz x nenergy
        """
        lz = self._zlayer(z)
        att = self.getcache("attenuationinfo")
        linatt = att["linatt"].loc[lz][energy]
        cor = att["linatt_cumulcor"].loc[lz][energy]
        if linatt.ndim != 0:
            linatt = linatt.values
            cor = cor.values
        if linatt.ndim == 2:
            z = z[:, np.newaxis]
        return z*linatt + cor

    def _transmission(self, zi, zj, cosaij, energy):
        """Transmission from depth zi to zj

        Args:
            zi(num|array): start depth of attenuation (nz)
            zj(num|array): end depth of attenuation (nz)
            cosaij(num|array): angle with surface normal (nz)
            energy(num|array): energies to be attenuation (nenergy)

        Returns:
            array: nz x nenergy
        """
        datt = self._cum_attenuation(
            zj, energy)-self._cum_attenuation(zi, energy)
        if datt.ndim == 2:
            if instance.isarray(cosaij):
                cosaij = cosaij[:, np.newaxis]
        # assert(sum(instance.asarray(-datt/cosaij)>0)==0)
        return np.exp(-datt/cosaij)

    def _cache_interactioninfo(self, energy, emin=None, emax=None, ninteractions=None, geomkwargs=None):
        # Pepare resulting lists
        _nlayers = self.nlayers+2
        _ninteractions = ninteractions+2
        probabilities = [None]*_ninteractions
        energyindex = [None]*_ninteractions

        # Source energies
        energy, func = instance.asarrayf(energy)
        nsource = len(energy)
        probabilities[0] = pd.DataFrame(
            columns=[xrayspectrum.RayleighLine(energy)])

        # Calculate interaction probabilities (ph/cm/srad)
        def getenergy(x, **kwargs):
            return list(listtools.flatten(interaction.energy(**kwargs) for interaction in x.columns))

        for i in range(ninteractions):
            energy = getenergy(probabilities[i], **geomkwargs)
            nenergy = len(energy)

            def f(x, energy=energy):
                return (np.abs(energy-x)).argmin()

            # Interaction probabilities (1/cm/srad):
            #  column -> interaction
            #  index -> [layer,source]
            probs = [None]*_nlayers
            probs[1:-1] = [pd.DataFrame.from_dict(dict(layer.xrayspectrum(
                energy, emin=emin, emax=emax).probabilities)) for layer in self]
            probs[0] = pd.DataFrame(index=range(nenergy))
            probs[-1] = probs[0]
            probs = pd.concat(probs, sort=True)
            probs.fillna(0., inplace=True)
            probs.index = pd.MultiIndex.from_product(
                [np.arange(_nlayers), range(nenergy)], names=["layer", "source"])

            probabilities[i+1] = probs
            energyindex[i+1] = f

        return {"probabilities": probabilities,
                "energyindex": energyindex,
                "getenergy": getenergy}

    def _genprobabilities(self, zi, i, energyi, interactionj):
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

        interactioninfo = self.getcache("interactioninfo")
        energyindex = interactioninfo["energyindex"][i](energyi)

        # Advanced indexing on MultiIndex: does not preserve order and repeats
        probs = interactioninfo["probabilities"][i].loc[(
            lz, energyindex), interactionj]
        if probs.ndim != 0:
            if lzarr:
                # apply order and repeats in lz (assume energyindex is a scalar)
                probs.index = probs.index.droplevel(1)
                probs = probs.loc[lz]
            probs = probs.values

        return probs

    def _gentransmission(self, zi, zj, cosaij, i, energyi, energyj, interactionj):
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
        probs = self._genprobabilities(zi, i, energyi, interactionj)
        T = self._transmission(zi, zj, cosaij, energyj)
        return probs*T

    def _gentransmission_saintegrated(self, zi, zj, i, energyi, energyj, interactionj):
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
        probs = self._genprobabilities(zi, i, energyi, interactionj)
        Aj = self._cum_attenuation(zj, energyj)
        Ai = self._cum_attenuation(zi, energyj)

        barri = instance.isarray(zi)
        barrj = instance.isarray(zj)
        if barri and barrj:
            probs = instance.asarray(probs)[:, np.newaxis]
            Ai = instance.asarray(Ai)[:, np.newaxis]
            Aj = instance.asarray(Aj)[np.newaxis, :]

        #func = lambda theta,phi: probs*np.exp(-(Aj-Ai)/np.cos(theta))*np.tan(theta)
        # return np.nquad(func,[(0,np.pi/2),(0,2*np.pi)])

        return (2*np.pi)*probs*scipy.special.exp1(Aj-Ai)

    def _primary_rates(self, selfabs=True):
        """Returns the ph generated per source line after 1 interaction (without efficiency term)
        """
        interactionindex = 1

        interactioninfo = self.getcache("interactioninfo")
        energy0 = interactioninfo["getenergy"](
            interactioninfo["probabilities"][0])

        nsource = len(energy0)
        interactions1 = interactioninfo["probabilities"][interactionindex].columns
        nlayers = self.nlayers
        nlines = len(interactions1)

        # Integrated attenuation over the sample thickness
        if selfabs:
            geomkwargs = self.geometry.xrayspectrumkwargs()
            energy1 = interactioninfo["getenergy"](
                interactioninfo["probabilities"][interactionindex], **geomkwargs)

            att = self.getcache("attenuationinfo")
            cosafirst = self.geometry.cosnormin
            cosalast = self.geometry.cosnormout
            mu0 = att["linatt"][energy0].values/cosafirst
            mu1 = att["linatt"][energy1].values/cosalast
            cor0 = att["linatt_cumulcor"][energy0].values/cosafirst
            cor1 = att["linatt_cumulcor"][energy1].values/cosalast

            chi = mu1[1:-1, :, np.newaxis] - mu0[1:-1, np.newaxis, :]
            chicor = cor1[1:-1, :, np.newaxis] - cor0[1:-1, np.newaxis, :]

            layerinfo = self.getcache("layerinfo")
            J2 = np.exp(chi*layerinfo["cumul_thickness"]
                        [1:, np.newaxis, np.newaxis])
            J2 -= np.exp(chi*layerinfo["cumul_thickness"]
                         [:-1, np.newaxis, np.newaxis])
            J2 /= chi
            J2 *= np.exp(chicor)
            if not self.geometry.reflection:
                J2 *= np.exp(-cor1[-1, np.newaxis, :, np.newaxis])

            # nlayers x nenergy1 x nenergy0 -> nlayers x nsource x nenergy1
            J2 = np.transpose(J2, [0, 2, 1])

            # nlayers x nenergy0 x nenergy1 -> nlayers x nsource x nlines (reduce scattering lines)
            interactions1exp = list(listtools.flatten(
                [interaction]*interaction.nenergy for interaction in interactions1))
            indC = np.asarray(
                [interaction == "Compton" for interaction in interactions1exp])
            indR = np.asarray(
                [interaction == "Rayleigh" for interaction in interactions1exp])
            indF = ~indC & ~indR

            indsource = range(nsource)
            J2 = np.concatenate((J2[..., indF],
                                 J2[:, indsource, indC][..., np.newaxis],
                                 J2[:, indsource, indR][..., np.newaxis]), axis=-1)
            interactions1 = interactions1.tolist()
            interactions1.insert(nlines, interactions1.pop(
                interactions1.index("Compton")))
            interactions1.insert(nlines, interactions1.pop(
                interactions1.index("Rayleigh")))
        else:
            # lim[x->0] (exp(x.d)-1)/x = d

            # nlayers x 1 x 1
            J2 = self.thickness[:, np.newaxis, np.newaxis]

        # cm -> cm.srad
        integratormult = self.geometry.solidangle/self.geometry.cosnormin
        J2 *= integratormult

        # Interaction probability: nlayers x nsource x nlines  (1/cm/srad)
        probs = interactioninfo["probabilities"][interactionindex].loc[(
            range(1, self.nlayers+1),), interactions1]
        probs = probs.values.reshape((nlayers, nsource, nlines))

        # Rate: fluoresence/scattering per incoming photon
        J2 = J2*probs  # ph/phsource

        # Sum over layers
        J2 = J2.sum(axis=0).T  # nlines x nsource

        # Sum over source lines for fluorescence
        J2 = dict(zip(interactions1, J2))

        return J2

    def _primary_rates_numerical(self):
        """Returns the ph generated per source line after 1 interaction (without efficiency term)
        """
        interactionindex = 1

        cosafirst = self.geometry.cosnormin
        cosalast = self.geometry.cosnormout
        integratormult = self.geometry.solidangle/cosafirst
        layerinfo = self.getcache("layerinfo")
        za = layerinfo["cumul_thickness"][0]
        zb = layerinfo["cumul_thickness"][-1]
        zfirst = layerinfo["cumul_thickness"][0]
        zlast = layerinfo["zexit"]
        geomkwargs = self.geometry.xrayspectrumkwargs()

        interactioninfo = self.getcache("interactioninfo")
        energy0 = interactioninfo["getenergy"](
            interactioninfo["probabilities"][0])
        interactions1 = interactioninfo["probabilities"][interactionindex].columns

        def numintegrate(path, za, zb):
            return scipy.integrate.quad(path, za, zb)[0]

        n = (zb-za)/min(self.thickness)*100

        def numintegratefast(path, za, zb):
            x = np.linspace(za, zb, n)
            y = path(x)
            return np.trapz(y, x=x)
            # return scipy.integrate.trapz(y, x=x)

        #import matplotlib.pyplot as plt
        J2 = {}
        for interaction1 in interactions1:
            energy1 = interaction1.energy(**geomkwargs)
            if isinstance(interaction1, xrayspectrum.FluoZLine):
                energy1 = [energy1]*len(energy0)
            gen = [numintegrate(
                lambda z1: self._transmission(zfirst, z1, cosafirst, en0) *
                self._gentransmission(z1, zlast, cosalast, interactionindex, en0, en1, interaction1), za, zb)
                for en0, en1 in zip(energy0, energy1)]
            # plt.figure()
            #x = np.linspace(za,zb,n)
            # plt.plot(x,path(x))
            # plt.show()
            J2[interaction1] = np.asarray(gen)*integratormult
        return J2

    def _secondary_interaction_numerical(self):
        """Returns the ph generated per source line after 2 interactions (without efficiency term)
        """
        interactionindex = 2

        cosafirst = self.geometry.cosnormin
        cosalast = self.geometry.cosnormout
        integratormult = self.geometry.solidangle/cosafirst
        layerinfo = self.getcache("layerinfo")
        za = layerinfo["cumul_thickness"][0]
        zb = layerinfo["cumul_thickness"][-1]
        zfirst = layerinfo["cumul_thickness"][0]
        zlast = layerinfo["zexit"]
        geomkwargs1 = self.geometry.xrayspectrumkwargs()
        geomkwargs2 = geomkwargs1

        interactioninfo = self.getcache("interactioninfo")
        energy0 = interactioninfo["getenergy"](
            interactioninfo["probabilities"][0])
        interactions1 = interactioninfo["probabilities"][1].columns
        interactions2 = interactioninfo["probabilities"][2].columns
        J3 = {}

        def path(z1, z2): return self._transmission(zfirst, z1, cosafirst, en0)[:, np.newaxis] *\
            self._gentransmission_saintegrated(z1, z2, interactionindex-1, en0, en1, interaction1) *\
            self._gentransmission(
                z2, zlast, cosalast, interactionindex, en1, en2, interaction2)[np.newaxis, :]

        def numintegrate(path, za, zb):
            return scipy.integrate.nquad(path, [(za, zb)]*2)[0]

        n = (zb-za)/min(self.thickness)*100

        def numintegratefast(path, za, zb):
            x1 = np.linspace(za, zb, n)
            x2 = np.linspace(za, zb, n)
            y = path(x1, x2)
            y = np.trapz(y, x=x1, axis=0)
            y = np.trapz(y, x=x2, axis=0)
            return y

        import matplotlib.pyplot as plt

        for interaction1 in interactions1:
            energy1 = interaction1.energy(**geomkwargs1)
            if isinstance(interaction1, xrayspectrum.FluoZLine):
                energy1 = [energy1]*len(energy0)

            for interaction2 in interactions2:
                energy2 = interaction2.energy(**geomkwargs2)
                if isinstance(interaction2, xrayspectrum.FluoZLine):
                    energy2 = [energy2]*len(energy1)

                for en0, en1, en2 in zip(energy0, energy1, energy2):
                    x1 = np.linspace(za, zb, n)
                    x2 = np.linspace(za, zb, n+1)

                    print(self._transmission(zfirst, x1, cosafirst, en0).shape)
                    print(self._gentransmission_saintegrated(
                        x1, x2, interactionindex-1, en0, en1, interaction1).shape)
                    print(self._gentransmission(x2, zlast, cosalast,
                                                interactionindex, en1, en2, interaction2).shape)

                    plt.figure()

                    img = path(x1, x2)
                    plt.imshow(img)
                    plt.show()

                gen = [numintegrate(path, za, zb)
                       for en0, en1, en2 in zip(energy0, energy1, energy2)]

                J3[interaction2] = np.asarray(gen)*integratormult

        return J3

    def addtofisx(self, setup, cfg):
        setup.setSample([layer.addtofisx(setup, cfg) for layer in self])
        self.geometry.addtofisx(setup, cfg)

    def addtopymca_matrix(self, setup, cfg, name, thickness=0.):
        anglein = self.geometry.anglein
        angleout = self.geometry.angleout
        scatteringangle = self.geometry.scatteringangle
        if name == "MULTILAYER":
            density = 0.
        else:
            v = cfg["materials"][name]
            density = v["Density"]
        cfg["attenuators"]["Matrix"] = [1, name, density,
                                        thickness, anglein, angleout, 0, scatteringangle]

    def loadfrompymca_matrix(self, setup, cfg):
        _, name, density, thickness, anglein, angleout, _, scatteringangle = cfg[
            "attenuators"]["Matrix"]
        self.geometry.anglein = anglein
        self.geometry.angleout = angleout
        return name, density, thickness

    def addtopymca_layer(self, setup, cfg, index, layer):
        name = setup.addtopymca_material(
            cfg, layer, defaultthickness=layer.thickness)
        l = "Layer{}".format(index)
        cfg["multilayer"][l] = [1, name, layer.density, layer.thickness]

    def loadfrompymca_layer(self, setup, cfg, index):
        l = "Layer{}".format(index)
        if l in cfg["multilayer"]:
            enabled, name, density, thickness = cfg["multilayer"][l]
            if enabled:
                material = setup.loadfrompymca_material(cfg, name, density)
                return (material, thickness)
            else:
                return tuple()
        else:
            return None

    def addtopymca_shells(self, setup, cfg, elements):
        emax = setup.emax_strict
        emin = setup.emin

        if "peaks" not in cfg:
            cfg["peaks"] = {}
        for e in elements:
            shells = e.pymcashellfactory(emin=emin, emax=emax)
            if shells:
                cfg["peaks"][str(e)] = shells

    def addtopymca(self, setup, cfg):
        if self.nlayers == 1:
            name = setup.addtopymca_material(
                cfg, self[0], defaultthickness=self[0].thickness)
            self.addtopymca_shells(setup, cfg, self[0].elements)
            self.addtopymca_matrix(
                setup, cfg, name, thickness=self[0].thickness)
        else:
            for index, layer in enumerate(self):
                self.addtopymca_layer(setup, cfg, index, layer)
                self.addtopymca_shells(setup, cfg, layer.elements)
            self.addtopymca_matrix(setup, cfg, 'MULTILAYER')
        self.geometry.addtopymca(setup, cfg)

    def loadfrompymca(self, setup, cfg):
        self.geometry.loadfrompymca(setup, cfg)
        name, density, thickness = self.loadfrompymca_matrix(setup, cfg)
        if name == "MULTILAYER":
            layer = tuple()
            index = 0
            layers = []
            while layer is not None:
                layer = self.loadfrompymca_layer(setup, cfg, index)
                index += 1
                if layer:
                    layers.append(layer)
            material, thickness = zip(*layers)
        else:
            material = [setup.loadfrompymca_material(cfg, name, density)]
            thickness = [thickness]

        self.layers = [Layer(material=mat, thickness=t, ml=self)
                       for mat, t in zip(material, thickness)]

    def _parse_fisx_result(self, gen):
        result = {}

        # Get fluorescence rates from fisx (add escape peaks)
        for group in gen:
            el = element.Element(group.split(' ')[0])
            for layer in gen[group]:
                for peak in gen[group][layer]:
                    line = xrayspectrum.FluoLine(peak.split(' ')[0])
                    line = xrayspectrum.FluoZLine(el, line)
                    rate = gen[group][layer][peak]["rate"]
                    if line in result:
                        result[line] += rate
                    else:
                        result[line] = rate

        # Correction for detector in transmission
        if not self.geometry.reflection:
            for line in result:
                energy = line.energy(**self.geometry.xrayspectrumkwargs())
                result[line] = result[line]*self.transmission(energy, out=True)

        return result

    def _dict_to_spectrum(self, gen, emin=0, emax=None, scattering=True):
        if not scattering:
            gen = {k: v for k, v in gen.items() if isinstance(
                k, xrayspectrum.FluoZLine)}

        if emax is None:
            allenergies = list(listtools.flatten(line.energy(
                **self.geometry.xrayspectrumkwargs()) for line in gen))
            emax = max(allenergies)

        spectrum = xrayspectrum.Spectrum()
        spectrum.update(gen)
        spectrum.xlim = [emin, emax]
        spectrum.density = None
        spectrum.title = str(self)
        spectrum.type = spectrum.TYPES.rate
        spectrum.geometry = self.geometry
        return spectrum

    def _print_fisx(self, fluo, details=False):
        if details:
            rowfmt = "{:>6}{:>8}{:>20}{:>10}{:>10}{:>20}{:>20}{:>20}{:>20}{:>20}"
            print(rowfmt.format("Layer", "Element", "MassFrac", "Line", "Energy",
                                "Rate", "Primary", "Multiplier(2)", "Multiplier(2+3)", "Efficiency"))
        else:
            rowfmt = "{:>6}{:>8}{:>20}{:>10}{:>10}{:>20}"
            print(rowfmt.format("Layer", "Element",
                                "MassFrac", "Line", "Energy", "Rate"))

        for key in sorted(fluo):
            ele = key.split(' ')[0]
            for layer in fluo[key]:
                lines = sorted(list(fluo[key][layer].keys()))
                for line in lines:
                    if line.endswith("esc"):
                        continue
                    # Mass fraction in this layer
                    w = fluo[key][layer][line]["massFraction"]
                    # energy of the line
                    energy = fluo[key][layer][line]["energy"]
                    # expected measured rate (everything except flux*time)
                    rate = fluo[key][layer][line]["rate"]
                    escaperate = sum(fluo[key][layer][line2]["rate"] for line2 in lines if line2.endswith(
                        "esc") and line2.startswith(line))
                    rate += escaperate
                    # primary photons (no attenuation and no detector considered)
                    primary = fluo[key][layer][line]["primary"]
                    # secondary photons (no attenuation and no detector considered)
                    secondary = fluo[key][layer][line]["secondary"]
                    # tertiary photons (no attenuation and no detector considered)
                    tertiary = fluo[key][layer][line].get("tertiary", 0.0)
                    # attenuation and detector
                    efficiency = fluo[key][layer][line].get("efficiency", 0.0)
                    # correction due to secondary excitation
                    enhancement2 = (primary + secondary) / primary
                    # correction due to tertiary excitation
                    enhancement3 = (primary + secondary + tertiary) / primary

                    if details:
                        print(rowfmt.format(layer, ele, w, line, energy, rate+escaperate,
                                            primary, enhancement2, enhancement3, efficiency))
                    else:
                        print(rowfmt.format(layer, ele, w,
                                            line, energy, rate+escaperate))

                    assert(np.isclose(
                        rate, (primary + secondary + tertiary)*efficiency))

    @staticmethod
    def _parse_weights(weights, nsource):
        if weights is None:
            weights = np.ones(nsource, dtype=float)/nsource
        else:
            weights = instance.asarray(weights)
            weights = weights/weights.sum(dtype=float)
        return weights

    def _rates_fisx(self, energy0, weights, ninteractions, emin=0, emax=None):
        energy0 = instance.asarray(energy0)
        nsource = len(energy0)
        if emax is None:
            emax = np.max(energy0)

        setup = fisx.XRF()

        # Add sample, detector and geometry
        self.addtofisx(setup, self.FISXCFG)

        # Add source
        setup.setBeam(energy0, weights=self._parse_weights(weights, nsource))

        # Peak groups
        groups = {}
        for layer in self:
            groups.update(layer.fisxgroups(emin=emin, emax=emax))

        def shellparse(shell):
            shell = str(shell)
            if shell.startswith("M"):
                shell = "M"
            return shell

        groups = {"{} {}".format(el, shellparse(shell))
                      for el, shells in groups.items() for shell in shells}

        # Get fluorescence
        secondary = 2*(ninteractions > 1)
        # 0: none, 1: intralayer, 2: interlayer
        gen = setup.getMultilayerFluorescence(
            groups, self.FISXCFG.FISXMATERIALS,
            secondary=secondary, useMassFractions=1)

        # self._print_fisx(gen)
        result = self._parse_fisx_result(gen)

        return result

    def _interactions_applyefficiency(self, gen, withdetectorattenuation=True):
        """Apply filter attenuation (source and detection) + detector efficiency
        """

        geom = self.geometry.xrayspectrumkwargs()

        lines = list(gen.keys())
        energysource = lines[lines.index("Rayleigh")].energy(**geom)
        energydet = [k.energy(**geom) for k in lines]
        ind = np.cumsum([listtools.length(en) for en in energydet])
        ind = np.insert(ind, 0, 0)
        ind = zip(ind[:-1], ind[1:])

        energydet = list(listtools.flatten(energydet))
        efficiency = self.geometry.efficiency(
            energysource, energydet, withdetectorattenuation=withdetectorattenuation)

        for k, (a, b) in zip(lines, ind):
            if a+1 == b:  # Fluorescence
                eff = efficiency[:, a]
            else:  # Scattering
                eff = np.diag(efficiency[:, a:b])
            gen[k] = gen[k]*eff

    @cache.withcache("layerinfo")
    def xrayspectrum(self, energy0, emin=0, emax=None, method="analytical",
                     ninteractions=1, weights=None, scattering=True,
                     withdetectorattenuation=True):

        if method == "fisx":
            if scattering:
                method = "analytical"
            # if not self.geometry.reflection and self.nlayers>1:
            #    method="analytical"
            if not self.geometry.reflection:
                method = "analytical"

            if ninteractions > 3:
                method = "analytical"

        if method == "analytical":
            if ninteractions >= 2:
                method = "numerical"

        if method == "fisx":
            gen = self._rates_fisx(
                energy0, weights, ninteractions, emin=emin, emax=emax)
        else:
            geomkwargs = self.geometry.xrayspectrumkwargs()
            with self.cachectx("interactioninfo", energy0, emin=emin, emax=emax,
                               ninteractions=ninteractions,
                               geomkwargs=geomkwargs):
                interactioninfo = self.getcache("interactioninfo")
                allenergies = interactioninfo["getenergy"](
                    interactioninfo["probabilities"][-2], **geomkwargs)
                with self.cachectx("attenuationinfo", allenergies):
                    # Primary interaction (with self-absorption)
                    if method == "numerical":
                        gen = self._primary_rates_numerical()
                    else:
                        gen = self._primary_rates()

                    # Secondary interaction (with self-absorption)
                    if ninteractions >= 2 and False:
                        for k, v in self._secondary_interaction_numerical().items():
                            if k in gen:
                                gen[k] += v
                            else:
                                gen[k] = v

            # Apply filter attenuation (source and detection) + detector efficiency
            self._interactions_applyefficiency(
                gen, withdetectorattenuation=withdetectorattenuation)

        spectrum = self._dict_to_spectrum(
            gen, emin=emin, emax=emax, scattering=scattering)

        if method != "fisx":
            spectrum.apply_weights(weights)
            # spectrum.sum_sources() # this is already done when fisx is used

        return spectrum

    def propagate(self, N, energy, interaction="transmission", forward=True):
        """Error propagation of a number of photons.

        Args:
            N(num|array): incomming number of photons with uncertainties
            energy(num|array): energies

        Returns:
            num|numpy.array
        """
        # Bernouilli processes: compounding is the same as multiplication
        #                       so we can multiply the probabilities
        if interaction == "transmission":
            probsuccess = self.transmission(energy)
        else:
            raise RuntimeError("{} not implemented yet".format(interaction))

        N, probsuccess = self.propagate_broadcast(N, probsuccess)

        if instance.isuscalar(N):
            process = noisepropagation.bernouilli(probsuccess)
            Nout = noisepropagation.compound(N, process, forward=forward)
        else:
            if forward:
                Nout = N*probsuccess
            else:
                Nout = N/probsuccess

        return Nout


factory = Multilayer.factory
registry = Multilayer.clsregistry

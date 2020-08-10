# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import scipy.integrate
import scipy.special
import collections
import fisx
import logging
from contextlib import contextmanager

from ..utils import instance
from ..utils import cache
from ..utils import listtools
from ..math import fit1d
from ..math.utils import weightedsum
from . import xrayspectrum
from ..simulation.classfactory import with_metaclass
from ..simulation import xrmc
from ..simulation import xmimsim
from ..math import noisepropagation
from . import pymca
from . import element
from ..materials import compoundfromdb
from ..materials import mixture
from ..materials import types
from ..utils.copyable import Copyable
from .utils import reshape_spectrum_lines
from ..io import localfs
from ..io import spe

logger = logging.getLogger(__name__)


class Layer(Copyable):
    def __init__(self, material=None, thickness=None, fixed=False, parent=None):
        """
        Args:
            material(compound|mixture|str): material composition
            thickness(num): thickness in cm
            fixed(bool): thickness and composition are fixed
            parent(Multilayer): part of this ensemble
        """
        if instance.isstring(material):
            ret = compoundfromdb.factory(material)
            if ret is None:
                raise RuntimeError("Invalid material {}".format(material))
            material = ret
        self.material = material
        self.thickness = thickness
        self.fixed = fixed
        self.parent = parent

    def __getstate__(self):
        return {
            "material": self.material,
            "thickness": self.thickness,
            "fixed": self.fixed,
        }

    def __setstate__(self, state):
        self.material = state["material"]
        self.thickness = state["thickness"]
        self.fixed = state["fixed"]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (
                self.material == other.material
                and self.thickness == other.thickness
                and self.fixed == other.fixed
            )
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __getattr__(self, attr):
        return getattr(self.material, attr)

    def __str__(self):
        return "{} um ({})".format(self.thickness * 1e4, self.material)

    @property
    def xraythicknessin(self):
        return self.thickness / self.parent.geometry.cosnormin

    @xraythicknessin.setter
    def xraythicknessin(self, value):
        self.thickness = value * self.parent.geometry.cosnormin

    @property
    def xraythicknessout(self):
        return self.thickness / self.parent.geometry.cosnormout

    @xraythicknessout.setter
    def xraythicknessout(self, value):
        self.thickness = value * self.parent.geometry.cosnormout

    def absorbance(self, energy, weights=None, out=False, **kwargs):
        kwargs.pop("decomposed", None)
        if out:
            thickness = self.xraythicknessout
        else:
            thickness = self.xraythicknessin
        return self.material.absorbance(energy, thickness, weights=weights, **kwargs)

    def addtofisx(self, setup, cfg):
        name = cfg.addtofisx_material(self.material)
        return [name, self.density, self.thickness]

    def fisxgroups(self, emin=0, emax=np.inf):
        return self.material.fisxgroups(emin=emin, emax=emax)

    def arealdensity(self):
        wfrac = self.material.elemental_massfractions()
        m = self.density * self.thickness
        return dict(zip(wfrac.keys(), np.asarray(list(wfrac.values())) * m))


class Multilayer(with_metaclass((Copyable, cache.Cache))):
    """
    Class representing a multilayer of compounds or mixtures
    """

    FISXCFG = pymca.FisxConfig()

    def __init__(
        self, material=None, thickness=None, fixed=False, geometry=None, name=None
    ):
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
            fixed = fixed * len(material)
        self.layers = [
            Layer(material=mat, thickness=t, fixed=f, parent=self)
            for mat, t, f in zip(material, thickness, fixed)
        ]
        if not name:
            name = "MULTILAYER"
        self.name = name
        super(Multilayer, self).__init__(force=True)

    def __getstate__(self):
        return {"layers": self.layers, "geometry": self.geometry}

    def __setstate__(self, state):
        self.layers = state["layers"]
        for layer in self.layers:
            layer.parent = self
        self.geometry = state["geometry"]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.layers == other.layers and self.geometry == other.geometry
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
        layers = "\n ".join(
            "Layer {}. {}".format(i, str(layer)) for i, layer in enumerate(self)
        )
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
        return {el: w / s for el, w in ret.items()}

    def change_elemental_massfraction(self, Z, wZ):
        # wZ * sum_li(w_il*rho_il*t_il) = sum_l(w_Zl*rho_Zl*t_zl)
        # a: layers that contain Z
        # b: layers that do not contain Z
        # wZ * sum_ai(w_ia*rho_a*t_a) + wZ * sum_bi(w_ib*rho_b*t_b) = sum_a(w_Za*rho_a*t_a) + sum_b(w_Zb*rho_b*t_b)
        #
        #   t_A = sum_a(t_a)
        #   t_B = sum_b(t_b)
        #   t_a = t_A*r_a
        #   t_b = t_B*r_b = t*r_b - t_A*r_b
        #   t_B = t - t_A
        #
        # denom = + wZ * sum_ai(w_ia*rho_a*r_a) - sum_a(w_Za*rho_a*r_a)
        #         - wZ * sum_bi(w_ib*rho_b*r_b) + sum_b(w_Zb*rho_b*r_b)
        # num = t * sum_b(w_Zb*rho_b*r_b) - wZ*t * sum_bi(w_ib*rho_b*r_b)
        # t_A = num/denom
        #
        #  w_Zb = 0
        #
        # num = t * wZ * sum_bi(w_ib*rho_b*r_b)
        # denom = sum_a(w_Za*rho_a*r_a) - wZ * [sum_bi(w_ib*rho_b*r_b) - sum_ai(w_ia*rho_a*r_a)]
        pass

    def elemental_molefractions(self):
        return self.mixlayers().elemental_molefractions()

    def elemental_equivalents(self):
        return self.mixlayers().elemental_equivalents()

    def mixlayers(self):
        n = len(self)
        if n == 0:
            return None
        elif n == 1:
            return self[0].material
        else:
            vfrac = self.thickness
            vfrac = vfrac / float(vfrac.sum())
            materials = [layer.material for layer in self]
            return mixture.Mixture(
                materials, vfrac, types.fraction.volume, name=self.name
            )

    def mass_att_coeff(self, energy):
        """Total mass attenuation coefficient

        Args:
            energy(num|array): keV
        Returns:
            array: nz x nenergy
        """
        return np.asarray(
            [instance.asarray(layer.mass_att_coeff(energy)) for layer in self]
        )

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
            return [
                layer.absorbance(energy, weights=weights, out=out, fine=fine)
                for layer in self
            ]
        else:
            return np.sum(
                [
                    layer.absorbance(energy, weights=weights, out=out, fine=fine)
                    for layer in self
                ],
                axis=0,
            )

    def transmission(
        self, energy, weights=None, out=False, fine=False, decomposed=False
    ):
        A = self.absorbance(
            energy, weights=weights, out=out, fine=fine, decomposed=decomposed
        )
        if decomposed:
            return A  # TODO: apply recursively
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
            return y / A[0]

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

    def _refinerhod(
        self, energy, absorbance, refinedattr, fixedattr, weights=None, **kwargs
    ):
        y = absorbance
        for layer in self.fixediter():
            y = y - layer.absorbance(energy)

        A = [layer.mass_att_coeff(energy) for layer in self.freeiter()]
        if weights is not None:
            A = [weightedsum(csi, weights=weights) for csi in A]

        params = self._refine_linear(A, y, **kwargs)
        for param, layer in zip(params, self.freeiter()):
            setattr(layer, refinedattr, param / getattr(layer, fixedattr))
            logger.info(
                'Refined {} of "{}": {}'.format(
                    refinedattr, layer, getattr(layer, refinedattr)
                )
            )

    def refinecomposition(
        self, energy, absorbance, weights=None, fixthickness=True, **kwargs
    ):
        y = absorbance
        for layer in self.fixediter():
            y = y - layer.absorbance(energy, weights=weights)

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
            w = w / s
            w = dict(zip(layer.parts.keys(), w))
            layer.change_fractions(w, "mass")

            if fixthickness:
                layer.density = s / layer.xraythicknessin
                logger.info(
                    'Refined density of "{}": {} g/cm^3'.format(layer, layer.density)
                )
            else:
                layer.xraythicknessin = s / layer.density
                logger.info(
                    'Refined thickness "{}": {} g/cm^3'.format(
                        layer, layer.xraythicknessin
                    )
                )

    def refinethickness(self, energy, absorbance, **kwargs):
        self._refinerhod(energy, absorbance, "xraythicknessin", "density", **kwargs)

    def refinedensity(self, energy, absorbance, **kwargs):
        self._refinerhod(energy, absorbance, "density", "xraythicknessin", **kwargs)

    def _cache_layerinfo(self):
        t = np.empty(self.nlayers + 1)
        np.cumsum(self.thickness, out=t[1:])
        t[0] = 0
        if self.geometry.reflection:
            zexit = 0.0
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
        linattout = np.empty((self.nlayers + 2, nenergies), dtype=linatt.dtype)
        linattout[1:-1, :] = linatt
        linattout[[0, -1], :] = 0  # outside sample (vacuum)

        # Cumulative linear attenuation coefficient (= linatt*z + correction)
        attall = (linatt * thickness).sum(axis=0)
        cor = np.empty((self.nlayers + 2, nenergies), dtype=attall.dtype)
        cor[0, :] = 0  # before sample (vacuum)
        cor[-1, :] = attall  # after sample

        for i in range(nenergies):
            tmp = np.subtract.outer(linatt[:, i], linatt[:, i])
            tmp *= thickness
            cor[1:-1, i] = np.triu(tmp).sum(axis=0)

        linattout = pd.DataFrame(
            linattout, columns=energy, index=range(self.nlayers + 2)
        )
        cor = pd.DataFrame(cor, columns=energy, index=range(self.nlayers + 2))

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
        return z * linatt + cor

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
        datt = self._cum_attenuation(zj, energy) - self._cum_attenuation(zi, energy)
        if datt.ndim == 2:
            if instance.isarray(cosaij):
                cosaij = cosaij[:, np.newaxis]
        # assert(sum(instance.asarray(-datt/cosaij)>0)==0)
        return np.exp(-datt / cosaij)

    def _cache_interactioninfo(
        self, energy, emin=None, emax=None, ninteractions=None, geomkwargs=None
    ):
        """
        Args:
            energy(array): nSource x nSourceLines
        """

        def getenergy(x, **kwargs):
            return list(listtools.flatten(line.energy(**kwargs) for line in x.columns))

        # probabilities: list of pandas dataframes (one for each interaction)
        #                which saves the interaction probability of a layer
        #                at a particular energy
        #   column: line as a result of an interaction
        #   index:  [layer_index, energy_index]
        #   value:  interaction probability (1/cm/srad)
        # energy_to_index: list of functions (one for each interaction)
        #                  to get the energy_index closest to an energy
        _nlayers = self.nlayers + 2
        _ninteractions = ninteractions + 2
        probabilities = [None] * _ninteractions
        energy_to_index = [None] * _ninteractions
        interactioninfo = {
            "probabilities": probabilities,
            "energy_to_index": energy_to_index,
            "getenergy": getenergy,
        }

        # Interaction 0 has no probabilities
        # this is the source, not the result of an interaction
        source = [xrayspectrum.RayleighLine(energy)]
        probabilities[0] = pd.DataFrame(columns=source)

        # Calculate interaction probabilities (ph/cm/srad)
        for i in range(ninteractions):
            # Line energies after previous interaction
            energyi = getenergy(probabilities[i], **geomkwargs)
            nenergyi = len(energyi)
            # Interaction probabilities of each energy with each layer
            probs = [None] * _nlayers
            probs[1:-1] = [
                pd.DataFrame.from_dict(
                    dict(
                        layer.xrayspectrum(energyi, emin=emin, emax=emax).probabilities
                    )
                )
                for layer in self
            ]
            probs[0] = pd.DataFrame(index=range(nenergyi))
            probs[-1] = probs[0]
            probs = pd.concat(probs, sort=True)
            probs.fillna(0.0, inplace=True)
            probs.index = pd.MultiIndex.from_product(
                [np.arange(_nlayers), range(nenergyi)],
                names=["layer_index", "energy_index"],
            )
            probabilities[i + 1] = probs
            # Get energy_index closest to an energy
            energy_to_index[i + 1] = lambda x: (np.abs(energyi - x)).argmin()

        return interactioninfo

    def _prob_interaction(self, zi, i, energyi, interactionj):
        """
        Probability of interaction at depth zi

        Args:
            zi(num|array): one or more depths
            i(num): interaction order (1, 2, ...)
            energyi(num): energy of photon that interacts
            interactionj(object|array): one of more interactions

        Returns:
            array:
        """
        lz = self._zlayer(zi)
        lzarr = instance.isarray(lz)
        if lzarr:
            lz = lz.tolist()
        interactioninfo = self.getcache("interactioninfo")
        energy_index = interactioninfo["energy_to_index"][i](energyi)
        # Advanced indexing on MultiIndex: does not preserve order and repeats
        probs = interactioninfo["probabilities"][i].loc[
            (lz, energy_index), interactionj
        ]
        if probs.ndim != 0:
            if lzarr:
                # apply order and repeats in lz
                probs.index = probs.index.droplevel(1)
                probs = probs.loc[lz]
            probs = probs.values
        return probs

    def _prob_interaction_transmission(
        self, zi, zj, cosaij, i, energyi, energyj, interactionj
    ):
        """
        Probability of interaction at depth zi and reaching zj
        under a particular angle

        Args:
            zi(num|array): start depth of attenuation
            zj(num|array): end depth of attenuation
            cosaij(num|array): angle with surface normal
            i(num): interaction order (1, 2, ...)
            energyi(num): energy of photon that interacts
            energyj(num|array): energy of interactionj
            interactionj(object|array):

        Returns:
            array:
        """
        probs = self._prob_interaction(zi, i, energyi, interactionj)
        T = self._transmission(zi, zj, cosaij, energyj)
        return probs * T

    def _prob_interaction_transmission_saintegrated(
        self, zi, zj, i, energyi, energyj, interactionj
    ):
        """
        Total probability of interaction at depth zi and reaching depth zj

        Args:
            zi(num|array): start depth of attenuation
            zj(num|array): end depth of attenuation
            i(num): interaction order (1, 2, ...)
            energyi(num): energy of photon that interacts
            energyj(num): energy of interactionj
            interactionj(): energies to be attenuation

        Returns:
            array:
        """
        probs = self._prob_interaction(zi, i, energyi, interactionj)
        Aj = self._cum_attenuation(zj, energyj)
        Ai = self._cum_attenuation(zi, energyj)
        barri = instance.isarray(zi)
        barrj = instance.isarray(zj)
        if barri and barrj:
            probs = instance.asarray(probs)[:, np.newaxis]
            Ai = instance.asarray(Ai)[:, np.newaxis]
            Aj = instance.asarray(Aj)[np.newaxis, :]
        # Integrate over solid angle of emission from zi to zj (hemisphere)
        # TODO: assume isotropic emission from zi for now
        # func = lambda theta,phi: probs*np.exp(-(Aj-Ai)/np.cos(theta))*np.tan(theta)
        # return np.nquad(func,[(0,np.pi/2),(0,2*np.pi)])
        return (2 * np.pi) * probs * scipy.special.exp1(Aj - Ai)

    def _primary_rates(self, selfabs=True):
        """
        Returns the ph generated per source line after 1 interaction (without efficiency term)

        returns:
            dict: line: rates (nSourceLines)
        """
        interactionindex = 1

        interactioninfo = self.getcache("interactioninfo")
        energy0 = interactioninfo["getenergy"](interactioninfo["probabilities"][0])

        nsource = len(energy0)
        interactions1 = interactioninfo["probabilities"][interactionindex].columns
        nlayers = self.nlayers
        nlines = len(interactions1)

        # Effective sample thickness (corrected for attenuation)
        if selfabs:
            geomkwargs = self.geometry.xrayspectrumkwargs()
            energy1 = interactioninfo["getenergy"](
                interactioninfo["probabilities"][interactionindex], **geomkwargs
            )

            att = self.getcache("attenuationinfo")
            cosafirst = self.geometry.cosnormin
            cosalast = self.geometry.cosnormout
            mu0 = att["linatt"][energy0].values / cosafirst
            mu1 = att["linatt"][energy1].values / cosalast
            cor0 = att["linatt_cumulcor"][energy0].values / cosafirst
            cor1 = att["linatt_cumulcor"][energy1].values / cosalast

            chi = mu1[1:-1, :, np.newaxis] - mu0[1:-1, np.newaxis, :]
            chicor = cor1[1:-1, :, np.newaxis] - cor0[1:-1, np.newaxis, :]

            layerinfo = self.getcache("layerinfo")
            J2 = np.exp(chi * layerinfo["cumul_thickness"][1:, np.newaxis, np.newaxis])
            J2 -= np.exp(
                chi * layerinfo["cumul_thickness"][:-1, np.newaxis, np.newaxis]
            )
            J2 /= chi
            J2 *= np.exp(chicor)
            if not self.geometry.reflection:
                J2 *= np.exp(-cor1[-1, np.newaxis, :, np.newaxis])

            # nlayers x nenergy1 x nenergy0 -> nlayers x nsource x nenergy1
            J2 = np.transpose(J2, [0, 2, 1])

            # nlayers x nenergy0 x nenergy1 -> nlayers x nsource x nlines (reduce scattering lines)
            interactions1exp = list(
                listtools.flatten(
                    [interaction] * interaction.nenergy for interaction in interactions1
                )
            )
            indC = np.asarray(
                [interaction == "Compton" for interaction in interactions1exp]
            )
            indR = np.asarray(
                [interaction == "Rayleigh" for interaction in interactions1exp]
            )
            indF = ~indC & ~indR

            indsource = range(nsource)
            J2 = np.concatenate(
                (
                    J2[..., indF],
                    J2[:, indsource, indC][..., np.newaxis],
                    J2[:, indsource, indR][..., np.newaxis],
                ),
                axis=-1,
            )
            interactions1 = interactions1.tolist()
            interactions1.append(interactions1.pop(interactions1.index("Compton")))
            interactions1.append(interactions1.pop(interactions1.index("Rayleigh")))
        else:
            # lim[chi->0] (exp(chi.thickness)-1)/chi = thickness

            # nlayers x 1 x 1
            J2 = self.thickness[:, np.newaxis, np.newaxis]

        # Multiply thickness with geometrical factor: cm -> cm.srad
        J2 *= self.geometry.solidangle / self.geometry.cosnormin

        # Interaction probability: nlayers x nsource x nlines  (1/cm/srad)
        probs = interactioninfo["probabilities"][interactionindex].loc[
            (range(1, self.nlayers + 1),), interactions1
        ]
        probs = probs.values.reshape((nlayers, nsource, nlines))

        # Rate: fluoresence/scattering per incoming photon
        J2 = J2 * probs  # ph/phsource

        # Sum over layers
        J2 = J2.sum(axis=0).T  # nlines x nsource
        return dict(zip(interactions1, J2))

    def _primary_rates_numerical(self):
        """Returns the ph generated per source line after 1 interaction (without efficiency term)
        """
        interactionindex = 1

        cosafirst = self.geometry.cosnormin
        cosalast = self.geometry.cosnormout
        integratormult = self.geometry.solidangle / cosafirst
        layerinfo = self.getcache("layerinfo")
        za = layerinfo["cumul_thickness"][0]
        zb = layerinfo["cumul_thickness"][-1]
        zfirst = layerinfo["cumul_thickness"][0]
        zlast = layerinfo["zexit"]
        geomkwargs = self.geometry.xrayspectrumkwargs()

        interactioninfo = self.getcache("interactioninfo")
        energy0 = interactioninfo["getenergy"](interactioninfo["probabilities"][0])
        interactions1 = interactioninfo["probabilities"][interactionindex].columns

        def numintegrate(path, za, zb):
            return scipy.integrate.quad(path, za, zb)[0]

        n = (zb - za) / min(self.thickness) * 100

        def numintegratefast(path, za, zb):
            x = np.linspace(za, zb, n)
            y = path(x)
            return np.trapz(y, x=x)
            # return scipy.integrate.trapz(y, x=x)

        # import matplotlib.pyplot as plt
        J2 = {}
        for interaction1 in interactions1:
            energy1 = interaction1.energy(**geomkwargs)
            if isinstance(interaction1, xrayspectrum.FluoZLine):
                energy1 = [energy1] * len(energy0)

            def pathgen(en0, en1):
                return lambda z1: self._transmission(
                    zfirst, z1, cosafirst, en0
                ) * self._prob_interaction_transmission(
                    z1, zlast, cosalast, interactionindex, en0, en1, interaction1
                )

            paths = [pathgen(en0, en1) for en0, en1 in zip(energy0, energy1)]
            rates = [numintegrate(path, za, zb) for path in paths]
            # if interaction1 == 'Compton':
            #    plt.figure()
            #    x = np.linspace(za, zb, n)
            #    for path in paths:
            #        plt.plot(x, path(x))
            #    plt.show()
            J2[interaction1] = np.asarray(rates) * integratormult
        return J2

    def _secondary_interaction_numerical(self):
        """Returns the ph generated per source line after 2 interactions (without efficiency term)
        """
        # TODO: not finished
        interactionindex = 2

        cosafirst = self.geometry.cosnormin
        cosalast = self.geometry.cosnormout
        integratormult = self.geometry.solidangle / cosafirst
        layerinfo = self.getcache("layerinfo")
        za = layerinfo["cumul_thickness"][0]
        zb = layerinfo["cumul_thickness"][-1]
        zfirst = layerinfo["cumul_thickness"][0]
        zlast = layerinfo["zexit"]
        geomkwargs1 = self.geometry.xrayspectrumkwargs()
        geomkwargs2 = geomkwargs1

        interactioninfo = self.getcache("interactioninfo")
        energy0 = interactioninfo["getenergy"](interactioninfo["probabilities"][0])
        interactions1 = interactioninfo["probabilities"][1].columns
        interactions2 = interactioninfo["probabilities"][2].columns
        J3 = {}

        def path(z1, z2):
            return (
                self._transmission(zfirst, z1, cosafirst, en0)[:, np.newaxis]
                * self._prob_interaction_transmission_saintegrated(
                    z1, z2, interactionindex - 1, en0, en1, interaction1
                )
                * self._prob_interaction_transmission(
                    z2, zlast, cosalast, interactionindex, en1, en2, interaction2
                )[np.newaxis, :]
            )

        def numintegrate(path, za, zb):
            return scipy.integrate.nquad(path, [(za, zb)] * 2)[0]

        n = (zb - za) / min(self.thickness) * 100

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
                energy1 = [energy1] * len(energy0)

            for interaction2 in interactions2:
                energy2 = interaction2.energy(**geomkwargs2)
                if isinstance(interaction2, xrayspectrum.FluoZLine):
                    energy2 = [energy2] * len(energy1)

                for en0, en1, en2 in zip(energy0, energy1, energy2):
                    x1 = np.linspace(za, zb, n)
                    x2 = np.linspace(za, zb, n + 1)

                    print(self._transmission(zfirst, x1, cosafirst, en0).shape)
                    print(
                        self._prob_interaction_transmission_saintegrated(
                            x1, x2, interactionindex - 1, en0, en1, interaction1
                        ).shape
                    )
                    print(
                        self._prob_interaction_transmission(
                            x2,
                            zlast,
                            cosalast,
                            interactionindex,
                            en1,
                            en2,
                            interaction2,
                        ).shape
                    )

                    plt.figure()

                    img = path(x1, x2)
                    plt.imshow(img)
                    plt.show()

                rates = [
                    numintegrate(path, za, zb)
                    for en0, en1, en2 in zip(energy0, energy1, energy2)
                ]

                J3[interaction2] = np.asarray(rates) * integratormult

        return J3

    def addtofisx(self, setup, cfg):
        setup.setSample([layer.addtofisx(setup, cfg) for layer in self])
        self.geometry.addtofisx(setup, cfg)

    def addtopymca_matrix(self, setup, cfg, name, thickness=0.0):
        anglein = self.geometry.anglein
        angleout = self.geometry.angleout
        scatteringangle = self.geometry.scatteringangle
        if name == "MULTILAYER":
            density = 0.0
        else:
            v = cfg["materials"][name]
            density = v["Density"]
        cfg["attenuators"]["Matrix"] = [
            1,
            name,
            density,
            thickness,
            anglein,
            angleout,
            0,
            scatteringangle,
        ]

    def loadfrompymca_matrix(self, setup, cfg):
        _, name, density, thickness, anglein, angleout, _, scatteringangle = cfg[
            "attenuators"
        ]["Matrix"]
        self.geometry.anglein = anglein
        self.geometry.angleout = angleout
        return name, density, thickness

    def addtopymca_layer(self, setup, cfg, index, layer):
        name = setup.addtopymca_material(cfg, layer, defaultthickness=layer.thickness)
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
                cfg, self[0], defaultthickness=self[0].thickness
            )
            self.addtopymca_shells(setup, cfg, self[0].elements)
            self.addtopymca_matrix(setup, cfg, name, thickness=self[0].thickness)
        else:
            for index, layer in enumerate(self):
                self.addtopymca_layer(setup, cfg, index, layer)
                self.addtopymca_shells(setup, cfg, layer.elements)
            self.addtopymca_matrix(setup, cfg, "MULTILAYER")
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

        self.layers = [
            Layer(material=mat, thickness=t, parent=self)
            for mat, t in zip(material, thickness)
        ]

    def _parse_fisx_result(self, fisxresult):
        """
        Args:
            fisxresult(dict): group:dict(layer:dict(line):dict)
        Returns:
            dict: line: rate
        """
        # Get fluorescence rates from fisx (add escape peaks)
        rates = {}
        for group, layers in fisxresult.items():
            el = element.Element(group.split(" ")[0])
            for layer, peaks in layers.items():
                for peak, peakinfo in peaks.items():
                    line = xrayspectrum.FluoLine(peak.split(" ")[0])
                    line = xrayspectrum.FluoZLine(el, line)
                    rate = peakinfo["rate"]
                    if line in rates:
                        rates[line] += rate
                    else:
                        rates[line] = rate
        # Correction for detector in transmission
        # TODO: correct?
        if not self.geometry.reflection:
            for line in rates:
                energy = line.energy(**self.geometry.xrayspectrumkwargs())
                result[line] *= self.transmission(energy, out=True)
        return rates

    def _rates_to_spectrum(self, rates, emin=0, emax=None, scattering=True):
        """
        Args:
            rates(dict): line: rate
            emin(Optional(num)):
            emax(Optional(num)):
            scattering(Optional(bool)):
        Returns:
            xrayspectrum.Spectrum
        """
        if not scattering:
            rates = {
                k: v
                for k, v in rates.items()
                if not isinstance(k, xrayspectrum.ScatteringLine)
            }
        if emax is None:
            emax = max(
                listtools.flatten(
                    line.energy(**self.geometry.xrayspectrumkwargs()) for line in rates
                )
            )
        return xrayspectrum.Spectrum(
            rates,
            xlim=[emin, emax],
            density=None,
            title=str(self),
            type=xrayspectrum.Spectrum.TYPES.rate,
            geometry=self.geometry,
        )

    def _print_fisx(self, fluo, details=False):
        """
        Args:
            fluo(dict): group:dict(layer:dict(line):dict)
        """
        if details:
            rowfmt = "{:>6}{:>8}{:>20}{:>10}{:>10}{:>20}{:>20}{:>20}{:>20}{:>20}"
            print(
                rowfmt.format(
                    "Layer",
                    "Element",
                    "MassFrac",
                    "Line",
                    "Energy",
                    "Rate",
                    "Primary",
                    "Multiplier(2)",
                    "Multiplier(2+3)",
                    "Efficiency",
                )
            )
        else:
            rowfmt = "{:>6}{:>8}{:>20}{:>10}{:>10}{:>20}"
            print(
                rowfmt.format("Layer", "Element", "MassFrac", "Line", "Energy", "Rate")
            )
        for key in sorted(fluo):
            ele = key.split(" ")[0]
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
                    escaperate = sum(
                        fluo[key][layer][line2]["rate"]
                        for line2 in lines
                        if line2.endswith("esc") and line2.startswith(line)
                    )
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
                        print(
                            rowfmt.format(
                                layer,
                                ele,
                                w,
                                line,
                                energy,
                                rate + escaperate,
                                primary,
                                enhancement2,
                                enhancement3,
                                efficiency,
                            )
                        )
                    else:
                        print(
                            rowfmt.format(
                                layer, ele, w, line, energy, rate + escaperate
                            )
                        )
                    assert np.isclose(
                        rate, (primary + secondary + tertiary) * efficiency
                    )

    def _rates_fisx(self, energy0, weights, ninteractions, emin=0, emax=None):
        """
        Args:
            energy0(array): nSource x nLines
            weights(array): nSource x nLines
            ninteractions(num):
            emin(Optional(num)):
            emax(Optional(num)):
        Returns:
            list(dict): line: rate (nSource)
        """
        # Add sample, detector and geometry
        setup = fisx.XRF()
        cfg = self.FISXCFG
        self.addtofisx(setup, cfg)

        def shellparse(shell):
            shell = str(shell)
            if not shell.startswith("K") and not shell.startswith("L"):
                # Only K and L splitting supported
                shell = shell[0]
            return shell

        # Get fluorescence
        secondary = 2 * (ninteractions > 1)
        # 0: none, 1: intralayer, 2: interlayer
        rates = []
        for energy0i, weightsi in zip(energy0, weights):
            # Peak groups
            groups = {}
            if emax is None:
                emaxi = np.max(energy0i)
            else:
                emaxi = emax
            for layer in self:
                groups.update(layer.fisxgroups(emin=emin, emax=emaxi))
            groups = {
                "{} {}".format(el, shellparse(shell))
                for el, shells in groups.items()
                for shell in shells
            }
            # Add source
            setup.setBeam(energy0i, weights=weightsi)
            # Calculate fluorescence
            fixresult = setup.getMultilayerFluorescence(
                groups, cfg.FISXMATERIALS, secondary=secondary, useMassFractions=1
            )
            # self._print_fisx(fixresult)
            rates.append(self._parse_fisx_result(fixresult))
        return rates

    def _rates_calc(
        self,
        method,
        energy0,
        weights,
        ninteractions,
        emin=0,
        emax=None,
        withdetectorresponse=True,
    ):
        """
        Args:
            energy0(array): nSource x nLines
            weights(array): nSource x nLines
            ninteractions(num):
            emin(Optional(num)):
            emax(Optional(num)):
        Returns:
            list(dict): line: rate (nSource)
        """
        geomkwargs = self.geometry.xrayspectrumkwargs()
        rates = []
        for energy0i, weightsi in zip(energy0, weights):
            if emax is None:
                emaxi = np.max(energy0i)
            else:
                emaxi = emax
            with self.cachectx(
                "interactioninfo",
                energy0i,
                emin=emin,
                emax=emaxi,
                ninteractions=ninteractions,
                geomkwargs=geomkwargs,
            ):
                interactioninfo = self.getcache("interactioninfo")
                allenergies = interactioninfo["getenergy"](
                    interactioninfo["probabilities"][-2], **geomkwargs
                )
                with self.cachectx("attenuationinfo", allenergies):
                    # Primary interaction (with self-absorption)
                    if method == "numerical":
                        ratesi = self._primary_rates_numerical()
                    else:
                        ratesi = self._primary_rates()
                    # Secondary interaction (with self-absorption)
                    if ninteractions >= 2 and False:  # TODO
                        for k, v in self._secondary_interaction_numerical().items():
                            if k in ratesi:
                                ratesi[k] += v
                            else:
                                ratesi[k] = v
            # Attenuation of source and detected X-rays
            self._attenuated_rates(ratesi, withdetectorattenuation=withdetectorresponse)
            # Apply source weights
            for k in ratesi:
                ratesi[k] = ratesi[k] * weightsi
            rates.append(ratesi)
        return rates

    def _attenuated_rates(self, rates, withdetectorattenuation=True):
        """
        Apply various attenuations: source filter, detector filter, detector

        Args:
            rates(dict): line: rate
        """
        # Flat list of lines
        lines = list(rates.keys())

        # Source and detected lines
        geom = self.geometry.xrayspectrumkwargs()
        energysource = lines[lines.index("Rayleigh")].energy(**geom)
        energydet = [k.energy(**geom) for k in lines]
        ind = np.cumsum([listtools.length(en) for en in energydet])
        ind = np.insert(ind, 0, 0)
        ind = zip(ind[:-1], ind[1:])

        # Efficiency (nSource x nLines): filter and detector attenuation
        energydet = list(listtools.flatten(energydet))
        efficiency = self.geometry.efficiency(
            energysource, energydet, withdetectorattenuation=withdetectorattenuation
        )
        for k, (a, b) in zip(lines, ind):
            if a + 1 == b:  # Fluorescence
                eff = efficiency[:, a]
            else:  # Scattering
                eff = np.diag(efficiency[:, a:b])
            rates[k] = rates[k] * eff

    @contextmanager
    def _xrmc_context(
        self,
        flux=1,
        time=1,
        convoluted=False,
        pulseproctime=0,
        source_distance=1000,
        beamsize=1e-4,
    ):
        with localfs.temp(remove=False) as path:
            path.mkdir()
            # Units: keV, cm, degrees and sec
            world = xrmc.XrmcWorldBuilder(
                str(path), atmosphere=self.geometry.atmosphere
            )
            world.define_source(flux=flux, distance=source_distance, beamsize=beamsize)
            # Make sure the sample is larger than the beam footprint
            detdistance = self.geometry.distance.to("cm").magnitude
            samplesize = min(beamsize * 1000, detdistance)
            samplesize = min(samplesize, source_distance)
            # All layers the same size and beam goes through the sample center
            nlayers = self.nlayers
            dxs = [samplesize] * nlayers
            dys = [samplesize] * nlayers
            oxs = [0] * nlayers
            oys = [0] * nlayers
            for layer, dx, dy, ox, oy in zip(self, dxs, dys, oxs, oys):
                world.sample.add_layer(
                    material=layer.material,
                    thickness=layer.thickness,
                    dhor=dx,
                    dvert=dy,
                    ohor=ox,
                    overt=oy,
                )
            world.sample.polar = self.geometry.anglenormin
            world.sample.azimuth = self.geometry.sample_azimuth
            # Add beam filters
            dx = samplesize
            dy = samplesize
            ox = 0
            oy = 0
            for layer in self.geometry.beamfilters():
                world.source.add_layer(
                    material=layer["material"],
                    thickness=layer["thickness"],
                    dhor=dx,
                    dvert=dy,
                    ohor=ox,
                    overt=oy,
                    surface=10,
                )
            # Add detector
            activearea = self.geometry.detector.activearea.to("cm**2").magnitude
            mcagain = self.geometry.detector.mcagain
            polar = self.geometry.scatteringangle
            azimuth = self.geometry.detector_azimuth
            if convoluted:
                response = {
                    "material": self.geometry.detector.material,
                    "thickness": self.geometry.detector.thickness,
                    "noise": self.geometry.detector.mcanoise * 0
                    + 1e-10,  # to obtain line spectrum
                    "fano": self.geometry.detector.mcafano * 0
                    + 1e-10,  # to obtain line spectrum
                    "pulseproctime": pulseproctime,
                }
            else:
                response = {}
            world.add_xrfdetector(
                distance=detdistance,
                activearea=activearea,
                polar=polar,
                azimuth=azimuth,
                hoffset=0,
                voffset=0,
                emin=0,
                emax=1,
                ebinsize=mcagain,
                forcedetect=True,
                multiplicity=10,
                time=time,
                response=response,
            )
            # Add detector filters
            dhor, dvert = world.detector.pixelsize
            for layer in self.geometry.detectorfilters(include_atmosphere=False):
                world.detector.add_layer(
                    material=layer["material"],
                    thickness=-layer["thickness"],
                    dhor=dhor,
                    dvert=dvert,
                )
            yield world

    def _sourceflux(self, energies, samplesourcedist, sampleflux=1e10):
        """Convert flux on sample to flux of the source

        Args:
            energies(array):
            samplesourcedist(num): in cm
            sampleflux(num)
        Returns:
            array: flux for each source line
        """
        atmosphere = self.geometry.atmosphere
        if atmosphere:
            sourcelineflux = sampleflux / atmosphere.transmission(
                energies, samplesourcedist
            )
        else:
            sourcelineflux = np.full_like(energies, sampleflux)
        sourcelineflux /= len(sourcelineflux)
        for layer in self.geometry.beamfilters(include_atmosphere=False):
            sourcelineflux = sourcelineflux / layer["material"].transmission(
                energies, layer["thickness"]
            )
        return sourcelineflux

    def _rates_xrmc(
        self,
        energy0,
        weights,
        ninteractions,
        emin=0,
        emax=None,
        withdetectorresponse=True,
    ):
        """
        Args:
            energy0(array): nSource x nLines
            weights(array): nSource x nLines
            ninteractions(num):
            emin(Optional(num)):
            emax(Optional(num)):
            withdetectorresponse(Optional(bool))

        Returns:
            list(dict): line: rate (nSource)
        """
        rates = []
        with self._xrmc_context(flux=1, time=1, convoluted=True) as world:
            for energy0i, weightsi in zip(energy0, weights):
                # Sample flux to source lines
                fluxi = self._sourceflux(energy0i, world.source.distance)
                flux = fluxi.sum()
                weightsi = fluxi / flux
                world.spectrum.lines = [
                    [en, 0, fl * w] for en, w, fl in zip(energy0i, weightsi, fluxi)
                ]
                # Detector energy range
                if emax is None:
                    emaxi = np.max(energy0i) + 1
                else:
                    emaxi = emax
                # world.detector.emin = self.geometry.detector.mcazero
                world.detector.emin = emin
                world.detector.emax = emaxi
                world.detector.ebinsize = self.geometry.detector.mcagain
                # Run simulation
                interactions = (0,) + (10000,) * ninteractions
                world.finalize(interactions=interactions)
                if not world.simulate():
                    raise RuntimeError("Simulation failed")
                # TODO: xrmc issue #49
                data, info = world.detector.result(convoluted=True)
                mca = data.sum(axis=tuple(range(data.ndim - 1)))
                import matplotlib.pyplot as plt

                plt.figure()
                plt.plot(mca, label="xmrc")
                plt.legend()
                # Extract lines
                mask = mca != 0
                mca = mca[mask] / flux
                xenergy = info["xenergy"][mask]
                # TODO: response already included
                # if withdetectorresponse:
                #    mca = mca * self.geometry.detector.attenuation(xenergy)
                linesi = [xrayspectrum.Line(en) for en in xenergy]
                ratesi = dict(zip(linesi, mca))
                rates.append(ratesi)
        return rates

    def _rates_xmimsim(
        self,
        energy0,
        weights,
        ninteractions,
        emin=0,
        emax=None,
        withdetectorresponse=True,
        source_distance=100,
        beamsize=1e-4,
        runxrmc=False,
    ):
        """
        Args:
            energy0(array): nSource x nLines
            weights(array): nSource x nLines
            ninteractions(num):
            emin(Optional(num)):
            emax(Optional(num)):
            withdetectorresponse(Optional(bool))

        Returns:
            list(dict), bool: line: rate (nSource), convoluted
        """
        rates = []
        with localfs.temp(remove=False) as path:
            path.mkdir()
            for energy0i, weightsi in zip(energy0, weights):
                # Sample flux to source lines
                fluxi = self._sourceflux(energy0i, source_distance)
                flux = fluxi.sum()
                weightsi = fluxi / flux
                sample = self._xmimsim_sample()
                # Run simulation
                ph = pymca.PymcaHandle(
                    energy=energy0i,
                    weights=weightsi,
                    emin=emin,
                    emax=emax,
                    ninteractions=ninteractions,
                    flux=flux,
                    time=1,
                    sample=sample,
                )
                xmimsim.run(
                    str(path),
                    pymcahandle=ph,
                    source_distance=source_distance,
                    beamsize=beamsize,
                    has_atmosphere=bool(self.geometry.atmosphere),
                    runxrmc=runxrmc,
                )
                if runxrmc:
                    # TODO: xrmc issue #49
                    data, info = xrmc.loadxrmcresult_xmimsim(str(path), convoluted=True)
                    mca = data.sum(axis=tuple(range(data.ndim - 1)))
                else:
                    mca, info = xmimsim.loadxmimsimresult(str(path), convoluted=False)
                import matplotlib.pyplot as plt

                plt.figure()
                plt.plot(mca, label="xmimsim")
                plt.legend()
                # Extract lines
                mask = mca != 0
                mca = mca[mask] / flux
                xenergy = info["xenergy"][mask]
                # TODO: response already included
                # if withdetectorresponse:
                #    mca = mca * self.geometry.detector.attenuation(xenergy)
                linesi = [xrayspectrum.Line(en) for en in xenergy]
                ratesi = dict(zip(linesi, mca))
                rates.append(ratesi)
        return rates

    def _xmimsim_sample(self):
        # Add atmosphere layer which is thick enough to include source and detector
        if self.geometry.atmosphere:
            atm_thickness = max(source_distance, self.geometry.distance) * 2
            lst = [(atmosphere, atm_thickness)] + [
                (layer.material, layer.thickness) for layer in self
            ]
            material, thickness = zip(*lst)
            return self.__class__(
                material=material,
                thickness=thickness,
                geometry=self.geometry,
                name=self.name,
            )
        else:
            return self

    def _assert_rate_parameters(
        self, method, ninteractions=1, scattering=True, withdetectorresponse=True
    ):
        """
        Modify the method based on requested features
        """
        if method == "xrmc":
            if not xrmc.installed():
                raise RuntimeError("'xrmc' is not installed")
            if not withdetectorresponse:
                raise RuntimeError("'xrmc' cannot disable detector response")
        elif method == "xmimsim":
            if not xmimsim.installed():
                raise RuntimeError("'xmimsim' is not installed")
            if not withdetectorresponse:
                raise RuntimeError("'xmimsim' cannot disable detector response")
        elif method == "fisx":
            if scattering:
                raise RuntimeError("'fisx' does not support scattering")
            if not withdetectorresponse:
                raise RuntimeError("'fisx' cannot disable detector response")
            if not self.geometry.reflection:
                raise RuntimeError("'fisx' does not support transmission geometry")
            if ninteractions > 3:
                raise RuntimeError(
                    "'fisx' does not support {} interactions".format(ninteractions)
                )
        elif method == "analytical":
            if ninteractions >= 2:
                raise RuntimeError(
                    "'analytical' does not support {} interactions".format(
                        ninteractions
                    )
                )
        return method

    @cache.withcache("layerinfo")
    def xrayspectrum(
        self,
        energy0,
        emin=0,
        emax=None,
        method="analytical",
        ninteractions=1,
        weights=None,
        scattering=True,
        withdetectorresponse=True,
        **kwargs
    ):
        """
        Spectrum of this sample measured under the associated gemetry

        Args:
            energy0(array): nLines or nSource x nLines
            emin:
            emax:
            method:
            ninteractions:
            weights(array): nLines or nSource x nLines
            scattering(bool): include scattering peaks
            withdetectorresponse(bool):
        Returns:
            Spectrum or list(Spectrum)
        """
        self._assert_rate_parameters(
            method,
            ninteractions=ninteractions,
            scattering=scattering,
            withdetectorresponse=withdetectorresponse,
        )
        # Calculate line rate dictionary for each source
        energy0, weights, singlespectrum, singleline = reshape_spectrum_lines(
            energy0, weights=weights
        )
        if method == "fisx":
            rates = self._rates_fisx(
                energy0, weights, ninteractions, emin=emin, emax=emax, **kwargs
            )
        elif method == "xrmc":
            rates = self._rates_xrmc(
                energy0,
                weights,
                ninteractions,
                emin=emin,
                emax=emax,
                withdetectorresponse=withdetectorresponse,
                **kwargs
            )
        elif method == "xmimsim":
            rates = self._rates_xmimsim(
                energy0,
                weights,
                ninteractions,
                emin=emin,
                emax=emax,
                withdetectorresponse=withdetectorresponse,
                **kwargs
            )
        else:
            rates = self._rates_calc(
                method,
                energy0,
                weights,
                ninteractions,
                emin=emin,
                emax=emax,
                withdetectorresponse=withdetectorresponse,
                **kwargs
            )
        # X-ray spectrum for each source
        spectra = [
            self._rates_to_spectrum(rdict, emin=emin, emax=emax, scattering=scattering)
            for rdict in rates
        ]
        if singlespectrum:
            return spectra[0]
        else:
            return spectra

    def convoluted_xrayspectrum(
        self,
        energy0,
        emin=0,
        emax=None,
        method="analytical",
        ninteractions=1,
        weights=None,
        scattering=True,
        escape=True,
        pileup=True,
        flux=1e9,
        time=1,
        **kwargs
    ):
        """
        Spectrum of this sample measured under the associated gemetry

        Args:
            energy0(array): nLines or nSource x nLines
            emin:
            emax:
            method:
            ninteractions:
            weights(array): nLines or nSource x nLines
            scattering(bool): include scattering peaks

        Returns:
            tuple or list(Spectrum)
        """
        self._assert_rate_parameters(
            method,
            ninteractions=ninteractions,
            scattering=scattering,
            withdetectorresponse=True,
        )
        if method == "xrmc":
            pass
        elif method == "xmimsim":
            pass
        else:
            result = self.xrayspectrum(
                energy0,
                emin=emin,
                emax=emax,
                method=method,
                ninteractions=ninteractions,
                weights=weights,
                scattering=scattering,
                **kwargs
            )
            kwargs = {"fluxtime": flux * time, "histogram": True}
            if isinstance(result, list):
                result = [s.sumspectrum(**kwargs) for s in result]
            else:
                result = result.sumspectrum(**kwargs)
        return result

    def propagate(self, N, energy, interaction="transmission", forward=True):
        """
        Error propagation of transmitted number of photons.

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
                Nout = N * probsuccess
            else:
                Nout = N / probsuccess
        return Nout


factory = Multilayer.factory
registry = Multilayer.clsregistry

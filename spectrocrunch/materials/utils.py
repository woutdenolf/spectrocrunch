# -*- coding: utf-8 -*-

import numpy as np
from . import element
from ..utils import instance


def elemental_crosssections(dictin, dictout=None, w=1):
    if dictout is None:
        dictout = {}
    for k, v in dictin.items():
        if isinstance(v["cs"], dict):
            elemental_crosssections(v["cs"], w=v["w"], dictout=dictout)
        else:
            if k in dictout:
                dictout[k] += w*v["w"]*v["cs"]
            else:
                dictout[k] = w*v["w"]*v["cs"]
    return dictout


def reshape_spectrum_lines(energy, weights=None):
    """
    Args:
        energy(num or array): source energies
                              shape = nSourceLines or nSource x nSourceLines
        weights(Optional(num or array): source line weights
                                        shape=nSourceLines or nSource x nSourceLines
    Returns:
        tuple: energy(nSource x nSourceLines), weights(nSource x nSourceLines), singlesource
    """
    if instance.isquantity(energy):
        energy = energy.to('keV').magnitude
    singlesource = True
    if instance.isarraynot0(energy):
        singlesource = instance.asarray(energy).ndim < 2
    energy = np.atleast_1d(energy).astype(float)
    if energy.ndim == 1:
        energy = energy[np.newaxis, :]
    else:
        energy = np.atleast_2d(energy)
    if weights is None:
        weights = np.ones_like(energy)
    else:
        weights = np.atleast_1d(weights).astype(float)
        if weights.ndim == 1:
            weights = weights[np.newaxis, :]
        else:
            weights = np.atleast_2d(weights)
    if energy.shape != weights.shape:
        raise ValueError('Energy and Weights do not have the same shape (nSource x nSourceLines)')
    weights = weights/weights.sum(axis=-1, dtype=float)[:, np.newaxis]
    return energy, weights, singlesource

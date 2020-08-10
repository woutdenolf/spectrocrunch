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
                dictout[k] += w * v["w"] * v["cs"]
            else:
                dictout[k] = w * v["w"] * v["cs"]
    return dictout


def reshape_spectrum_lines(energy, weights=None, normalize=True, **others):
    """
    Args:
        energy(num or array): source energies
                              shape = [nSource x] nSourceLines
        weights(Optional(num or array): source line weights
                                        shape= [nSource x] nSourceLines
        **other: others to be reshaped
    Returns:
        tuple: energy(nSource x nSourceLines), weights(nSource x nSourceLines), *others, singlesource, singleline
    """
    # Reshape energy
    if instance.isquantity(energy):
        energy = energy.to("keV").magnitude
    shape = np.asarray(energy).shape
    if shape:
        singlesource = len(shape) == 1
        singleline = False
    else:
        singlesource = singleline = True
    energy = np.atleast_1d(energy).astype(float)
    if energy.ndim == 1:
        energy = energy[np.newaxis, :]
    else:
        energy = np.atleast_2d(energy)
    # Reshape and normalize weights
    if weights is None:
        weights = np.ones_like(energy) / float(energy.shape[1])
    else:
        weights = np.atleast_1d(weights).astype(float)
        if weights.ndim == 1:
            weights = weights[np.newaxis, :]
        else:
            weights = np.atleast_2d(weights)
        if normalize:
            weights = weights / weights.sum(axis=-1, dtype=float)[:, np.newaxis]
    if energy.shape != weights.shape:
        raise ValueError(
            "Energy and Weights do not have the same shape (nSource x nSourceLines)"
        )
    # Reshape others
    result = [energy, weights]
    for name, values in others.items():
        values = np.atleast_1d(values).astype(float)
        if values.ndim == 1:
            values = values[np.newaxis, :]
        else:
            values = np.atleast_2d(values)
        if energy.shape != values.shape:
            raise ValueError(
                "Energy and {} do not have the same shape (nSource x nSourceLines)".format(
                    name
                )
            )
    result += [singlesource, singleline]
    return result

# -*- coding: utf-8 -*-

from ..patch import xraylib
from ..utils import instance
import numpy as np
import pandas as pd


def eval(method, Z, E, applypost=True, dataframe=False):
    """
    Args:
        Z(array or num): atomic numbers
        E(array or num or Quantity): energy in keV
        applypost(Optional(bool)): shrink output dimensions to
                                   match input dimensions
        dataframe(Optional(bool)):
    Returns:
        array or num or Dataframe: shape = len(Z) x len(E), unless applypost is True
    """
    is_arr_Z = instance.isarray(Z)
    is_arr_E = instance.isarray(E)
    if instance.isquantity(E):
        E = E.to('keV').magnitude
    if is_arr_Z or is_arr_E:
        Z = np.atleast_1d(Z)
        E = np.atleast_1d(E).astype(float)
        shapeZ = Z.shape
        shapeE = E.shape
        Z = Z.flatten()
        E = E.flatten()
        if xraylib.xraylib_np:
            method = getattr(xraylib.xraylib_np, method)
            result = method(Z, E)
        else:
            method = getattr(xraylib, method)
            shape = Z.size, E.size
            result = np.empty(shape, dtype=float)
            for i, Zi in enumerate(Z):
                for j, Ej in enumerate(E):
                    result[i, j] = method(Zi, Ej)
        if is_arr_Z and is_arr_E:
            def postfunc(x):
                return x
        elif is_arr_Z:
            def postfunc(x):
                return x[:, 0].reshape(shapeZ)
        else:
            def postfunc(x):
                return x[0, :].reshape(shapeE)
    else:
        method = getattr(xraylib, method)
        result = method(Z, float(E))
        def postfunc(x):
            return x
    if dataframe:
        return pd.DataFrame(result, index=Z, columns=E)
    else:
        if applypost:
            return postfunc(result)
        else:
            return result, postfunc

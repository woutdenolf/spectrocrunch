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
        E = E.to("keV").magnitude
    if is_arr_Z or is_arr_E:
        is_array_Z = instance.isarraynot0(Z)
        is_array_E = instance.isarraynot0(E)
        Z = np.atleast_1d(Z)
        E = np.atleast_1d(E).astype(float)
        if is_array_Z:
            shapeZ = Z.shape
        else:
            shapeZ = tuple()
        if is_array_E:
            shapeE = E.shape
        else:
            shapeE = tuple()
        shape_result = shapeZ + shapeE

        def postfunc(x):
            return x.reshape(shape_result)

        Z = Z.flatten()
        E = E.flatten()
        result = None
        if xraylib.xraylib_np:
            cmethod = getattr(xraylib.xraylib_np, method)
            try:
                result = cmethod(Z, E)
            except ValueError:
                pass
        if result is None:
            cmethod = getattr(xraylib, method)
            shape = Z.size, E.size
            result = np.zeros(shape, dtype=float)
            for i, Zi in enumerate(Z):
                for j, Ej in enumerate(E):
                    try:
                        result[i, j] = cmethod(Zi, Ej)
                    except ValueError:
                        pass
    else:
        cmethod = getattr(xraylib, method)
        try:
            result = cmethod(Z, float(E))
        except ValueError:
            result = 0.0

        def postfunc(x):
            return x

    if dataframe:
        return pd.DataFrame(result, index=Z, columns=E)
    else:
        if applypost:
            return postfunc(result)
        else:
            return result, postfunc

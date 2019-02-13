# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
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

from ..patch import xraylib
from ..utils import instance
import numpy as np


def identity(x):
    return x


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
        Z = np.atleast_1d(Z).astype(int)
        E = np.atleast_1d(E).astype(float)
        if xraylib.xraylib_np:
            method = getattr(xraylib.xraylib_np, method)
            result = method(Z, E)
        else:
            method = getattr(xraylib, method)
            s = len(Z), len(E)
            result = np.empty(s, dtype=float)
            for i, Zi in enumerate(Z):
                for j, Ej in enumerate(E):
                    result[i, j] = method(Zi, Ej)
        if is_arr_Z and is_arr_E:
            postfunc = identity
        elif is_arr_Z:
            def postfunc(x): return x[:, 0]
        else:
            def postfunc(x): return x[0, :]
    else:
        method = getattr(xraylib, method)
        result = method(Z, E)
        postfunc = identity
    if dataframe:
        return pd.DataFrame(result, index=Z, columns=E)
    else:
        if applypost:
            return postfunc(result)
        else:
            return result, postfunc

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

from . import compound
from . import types
from ..utils import instance

import warnings
try:
    import xraylib
    registry = xraylib.GetCompoundDataNISTList()
except ImportError:
    xraylib = None
    warnings.warn("xraylib is not installed", ImportWarning)
    registry = []


class CompoundFromNist(compound.Compound):
    """Interface to a compound defined by a list of elements
    """

    def __init__(self,nistname,name=None):
        """
        Args:
            nistname(str): NIST name
            name(Optional[str]): compound name
        """

        data = xraylib.GetCompoundDataNISTByName(nistname)
        if name is None:
            name = data["name"]
        super(CompoundFromNist,self).__init__(data["Elements"],data["massFractions"],
                                        types.fraction.mass,data["density"],name=name)

def factory(name):
    return CompoundFromNist(name)

def search(name):
    name = name.lower()
    ret = [k for k in registry if name in k.lower()]
    if len(ret)>1:
        ret2 = [k for k in registry if name == k.lower()]
        if ret2:
            ret  = ret2
    if ret:
        return ret[0]
    else:
        return None

def compoundfromnist(name):
    return factory(name)
    

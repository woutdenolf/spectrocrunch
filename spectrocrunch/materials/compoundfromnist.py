# -*- coding: utf-8 -*-

from . import compound
from . import compoundsearch
from . import types
from ..utils import instance

import warnings
try:
    import xraylib
    def getnames():
        return list(xraylib.GetCompoundDataNISTList())
except ImportError:
    xraylib = None
    warnings.warn("xraylib is not installed", ImportWarning)
    def getnames():
        return []


class CompoundFromNist(compound.Compound):
    """Interface to a compound defined by a list of elements
    """

    def __init__(self, nistname, name=None):
        """
        Args:
            nistname(str): NIST name
            name(Optional[str]): compound name
        """

        data = xraylib.GetCompoundDataNISTByName(nistname)
        if name is None:
            name = data["name"]
        super(CompoundFromNist, self).__init__(data["Elements"], data["massFractions"],
                                               types.fraction.mass, data["density"], name=name)


def factory(name):
    return CompoundFromNist(name)


def search(name):
    return compoundsearch.search(getnames(), name)


def compoundfromnist(name):
    return factory(name)

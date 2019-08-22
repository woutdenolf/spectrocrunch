# -*- coding: utf-8 -*-

from . import base
from ..utils.classfactory import with_metaclass
from ..materials import multilayer
from ..geometries.base import FlatSample
from ..utils import instance


class XrayOptics(with_metaclass(base.Optics)):
    pass


class KB(XrayOptics):
    aliases = ["kb"]


class ZonePlate(XrayOptics):
    aliases = ["zp"]


class Filter(XrayOptics):
    aliases = ["filter"]

    def __init__(self, material=None, thickness=None):
        geom = FlatSample(anglein=90)
        self.multilayer = multilayer.Multilayer(material=material,
                                                thickness=thickness,
                                                geometry=geom)
        super(Filter, self).__init__(uselut=False)

    def __getstate__(self):
        state = super(Filter, self).__getstate__()
        state['multilayer'] = self.multilayer
        return state

    def __setstate__(self, state):
        self.multilayer = state['multilayer']
        super(Filter, self).__setstate__(state)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if not (super(Filter, self).__eq__(other)):
                return False
            return self.multilayer == other.multilayer
        else:
            return False

    def __str__(self):
        return "Filter {}".format(str(self.multilayer))

    def transmission(self, energy):
        return self.multilayer.transmission(energy)


factory = XrayOptics.factory
registry = XrayOptics.clsregistry

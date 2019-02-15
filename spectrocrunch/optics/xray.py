# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
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

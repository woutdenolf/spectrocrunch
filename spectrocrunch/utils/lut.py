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

from scipy import interpolate
import numpy as np

from . import instance
from . import units
from . import listtools


class LUT(object):

    def __init__(self,default=None,kind="linear"):
        self.clear(default=default)
        self.kind = kind

    def __getstate__(self):
        return {'kind': self.kind,
                '_x': self.x,
                '_y': self.y,
                '_default': self._default}

    def __setstate__(self, state):
        self.kind = state['kind']
        self.clear(default=state['_default'])
        self.x = state['_x']
        self.y = state['_y']

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if (self.x is None) ^ (other.x is None):
                return False
            if (self.y is None) ^ (other.y is None):
                return False
            if self.x is not None:
                if not (self.x == other.x).all():
                    return False
            if self.y is not None:
                if not (self.y == other.y).all():
                    return False
            return self.kind == other.kind and \
                   self._default == other._default
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        s = '\n '.join("{:~}: {:~}".format(*xy)
                       for xy in zip(self.x, self.y))
        if s:
            return "Lookup table:\n {}".format(s)
        else:
            return "Lookup table: {}".format(self(None))
            
    def clear(self, default=None):
        self.x = None
        self.y = None
        self._func = lambda x: default
        self._default = default
    
    def isempty(self):
        return self.x is None

    @property
    def xunits(self):
        if self.x is None:
            return units.ureg.dimensionless
        else:
            return self.x.units
    
    @property
    def yunits(self):
        if self.y is None:
            return units.ureg.dimensionless
        else:
            return self.y.units

    def zip(self, xunits, yunits):
        x, y = self.x, self.y
        if x is None:
            return zip([], [])
        if xunits:
            x = x.to(xunits)
        else:
            x = x.magnitude
        if yunits:
            y = y.to(yunits)
        else:
            y = y.magnitude
        return zip(x, y)

    def __call__(self, x):
        x, func = units.asqarrayf(x)
        x = x.to(self.xunits).magnitude
        y = self._func(x)
        y = units.Quantity(y, units=self.yunits)
        return func(y)

    def add(self, x, y):
        x = units.asqarray(x)
        y = units.asqarray(y)
        if self.x is None:
            self.x = x
            self.y = y
            x = x.magnitude
            y = y.magnitude
        else:
            x = x.to(self.xunits).magnitude
            y = y.to(self.yunits).magnitude
            x = np.append(self.x.magnitude, x)
            y = np.append(self.y.magnitude, y)
        x, y = listtools.sort2lists(x, y)
        self.x = units.Quantity(x, self.xunits)
        self.y = units.Quantity(y, self.yunits)

        if len(x)==1:
            y = y[0]
            self._func = lambda x: y
        else:
            # units get lost
            self._func = interpolate.interp1d(x,y,bounds_error=False,fill_value=(y[0],y[-1]),kind=self.kind)

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
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from copy import copy

from . import instance
from . import units
from . import listtools


class LUT(object):
    """Lookup table with sorted/unique keys and keys/values with units
    """

    def __init__(self, default=np.nan, kind="linear", x=None, y=None):
        self.clear(default=default)
        self.kind = kind
        if x is not None:
            if y is None:
                try:
                    y = [self.default.magnitude]*len(x)
                except TypeError:
                    y = self.default.magnitude
                yunits = self.default.units
                y = units.Quantity(y, yunits)
            self.add(x, y)

    def __getstate__(self):
        return {'kind': self.kind,
                '_x': self.x,
                '_y': self.y,
                'default': self.default}

    def __setstate__(self, state):
        self.kind = state['kind']
        self.clear(default=state['default'])
        self.x = state['_x']
        self.y = state['_y']

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if (self.x is None) ^ (other.x is None):
                return False
            if (self.y is None) ^ (other.y is None):
                return False
            if self.x is not None:
                if not np.array_equal(self.x, other.x):
                    return False
            if self.y is not None:
                if not np.array_equal(self.y, other.y):
                    return False
            if np.isnan(self.default) ^ np.isnan(other.default):
                return False
            if not np.isnan(self.default):
                if self.default != other.default:
                    return False
            return self.kind == other.kind
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        s = '\n '.join(["{:~}: {:~}".format(*xy) for xy in self])
        if s:
            return "Lookup table:\n {}".format(s)
        else:
            return "Lookup table: {}".format(self(None))

    def clear(self, default=None):
        self.x = None
        self.y = None
        self._func = None
        self.default = units.Quantity(default, units.ureg.dimensionless)

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

    def __iter__(self):
        for xy in self.zip(self.xunits, self.yunits):
            yield xy

    def items(self):
        return iter(self)

    def __len__(self):
        try:
            return len(self.x)
        except TypeError:
            return 1

    def __call__(self, x):
        x, func = units.asqarrayf(x)
        if self.isempty():
            y = [self.default.magnitude]*len(x)
            yunits = self.default.units
        else:
            x = x.to(self.xunits).magnitude
            y = self._func(x)
            yunits = self.yunits
        y = units.Quantity(y, units=yunits)
        return func(y)

    def _new_data(self, args, add=True):
        nargs = len(args)
        if nargs == 1:
            xy = args[0]
            if isinstance(xy, self.__class__):
                x, y = xy.x, xy.y
            elif isinstance(xy, (list, tuple)):
                x, y = xy
            else:
                raise ValueError
        elif nargs == 2:
            x, y = args
        else:
            raise TypeError('Expects one or two positional arguments')
        x = units.asqarray(x)
        y = units.asqarray(y)
        if self.isempty():
            self.x = x
            self.y = y
            x = x.magnitude
            y = y.magnitude
        else:
            x = x.to(self.xunits).magnitude
            y = y.to(self.yunits).magnitude
            x = np.append(x, self.x.magnitude)
            y = np.append(y, self.y.magnitude)
        x, y = listtools.unique2lists(x, y, add=add)
        x, y = listtools.sort2lists(x, y)
        self.x = units.Quantity(x, self.xunits)
        self.y = units.Quantity(y, self.yunits)
        if len(x) == 1:
            y = [y[0]]
            self._func = lambda x: y*len(x)
        else:
            self._func = interpolate.interp1d(x, y, bounds_error=False,
                                              fill_value=(y[0], y[-1]),
                                              kind=self.kind)

    def add(self, *args):
        self._new_data(args, add=True)

    def replace(self, *args):
        self._new_data(args, add=False)

    def __add__(self, xy):
        o = copy(self)
        o.add(xy)
        return o

    def __iadd__(self, xy):
        self.add(xy)
        return self

    def plot(self, **kwargs):
        if self.isempty():
            x = self.x.magnitude
            y = self.y.magnitude
            plt.plot(x, y, **kwargs)
            plt.xlabel('{:~}'.format(self.xunits))
            plt.ylabel('{:~}'.format(self.yunits))
        else:
            plt.axhline(y=self.default.magnitude)
            plt.xlabel('{:~}'.format(self.xunits))
            plt.ylabel('{:~}'.format(self.yunits))

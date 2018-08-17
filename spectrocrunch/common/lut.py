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

from operator import itemgetter
from scipy import interpolate
import numpy as np

from . import instance
from . import units
from ..patch.pint import ureg

class LUT(object):

    def __init__(self,default=None,interpolation="linear"):
        self.clear(default=default)
        self.kind = interpolation

    def __str__(self):
        s = '\n '.join("{}: {}".format(k,v) for k,v in self.table())
        if s:
            return "Lookup table:\n {}".format(s)
        else:
            return "Lookup table: {}".format(self(None))
            
    def clear(self,default=None):
        self._tbl = {}
        self._func = lambda x: default
        self._xunits = ureg.dimensionless
        self._yunits = ureg.dimensionless
        
    def __call__(self,x):
        x = units.asqarray(x).to(self._xunits).magnitude
        y = self._func(x)
        if self._yunits != ureg.dimensionless:
            y = units.Quantity(y,units=self._yunits)
        return y
        
    def table(self):
        return sorted(self._tbl.items(), key=itemgetter(0))
    
    def isempty(self):
        return not bool(self._tbl)
    
    def add(self,x,y):
        x = instance.asarray(x)
        y = instance.asarray(y)
        self._tbl.update(dict(zip(x,y)))

        (x,y) = zip(*self.table())
        x = units.asqarray(x)
        y = units.asqarray(y)
        self._xunits = x.units
        self._yunits = y.units
            
        if len(x)==1:
            y = y[0].magnitude
            self._func = lambda x: y
        else:
            # units get lost
            self._func = interpolate.interp1d(x,y,bounds_error=False,fill_value=(y[0],y[-1]),kind=self.kind)


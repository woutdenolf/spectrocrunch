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

from scipy import interpolate

from operator import itemgetter

import numpy as np

from ..common import instance

class KB(object):

    def __init__(self):
        self._tbl = {}
        self._transmission = lambda energy: 1
    
    def __str__(self):
        s = '\n '.join("{} keV: {} %".format(k,v*100) for k,v in sorted(self._tbl.items()))
        return "KB transmission:\n {}".format(s)
    
    def transmission(self,energy):
        return self._transmission(energy)
    
    def set_transmission(self,energy,transmission):
        if instance.isnumber(energy):
            self._tbl[energy] = transmission
        else:
            self._tbl.update(dict(zip(energy,transmission)))
 
        (x,y) = zip(*sorted(self._tbl.items(), key=itemgetter(0)))

        if len(x)==1:
            self._transmission = lambda energy: y[0]
        else:
            self._transmission = interpolate.interp1d(list(x),list(y),bounds_error=False,fill_value=(x[0],x[-1]))


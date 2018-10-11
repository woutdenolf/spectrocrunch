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

import numpy as np
import sympy
from sympy.utilities.lambdify import lambdify, implemented_function

from ..utils import instance

def eval(expr,subs):
    for x,v in subs.items():
        expr = expr.subs(x,v)
    return expr.evalf()


class clip(sympy.Function):

    def _eval_evalf(self, prec):
        return np.clip(*self.args) 

    def inverse(self, argindex=1):
        return iclip


class iclip(sympy.Function):

    def _eval_evalf(self, prec):
        x,cmin,cmax = self.args
        y,func = instance.asarrayf(x)
        y[y<cmin] = np.nan
        y[y>cmax] = np.nan
        return func(y)

    def inverse(self, argindex=1):
        return clip




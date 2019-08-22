# -*- coding: utf-8 -*-

import numpy as np
import sympy
from sympy.utilities.lambdify import lambdify, implemented_function

from ..utils import instance


def eval(expr, subs):
    for x, v in subs.items():
        expr = expr.subs(x, v)
    return expr.evalf()


class clip(sympy.Function):

    def _eval_evalf(self, prec):
        return np.clip(*self.args)

    def inverse(self, argindex=1):
        return iclip


class iclip(sympy.Function):

    def _eval_evalf(self, prec):
        x, cmin, cmax = self.args
        y, func = instance.asarrayf(x)
        y[y < cmin] = np.nan
        y[y > cmax] = np.nan
        return func(y)

    def inverse(self, argindex=1):
        return clip

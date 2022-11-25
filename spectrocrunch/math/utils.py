# -*- coding: utf-8 -*-

from __future__ import division

import numpy as np
import math
import fractions
import functools
from ..utils import instance


def logscale(img):
    ret = -np.log(img / np.nanmax(img))
    ret /= np.nanmax(ret)
    return 1 - ret


def sig_to_ndigits(x, sig):
    return sig - int(math.floor(math.log10(abs(x)))) - 1


round_ndigits = np.round


def ceil_ndigits(x, n):
    m = 10**n
    return math.ceil(x * m) / m


def floor_ndigits(x, n):
    m = 10**n
    return math.floor(x * m) / m


def round_sig(x, sig):
    return round_ndigits(x, sig_to_ndigits(x, sig))


def ceil_sig(x, sig):
    return ceil_ndigits(x, sig_to_ndigits(x, sig))


def floor_sig(x, sig):
    return floor_ndigits(x, sig_to_ndigits(x, sig))


def floatformat(x, sig):
    n = max(sig_to_ndigits(x, sig), 0)
    y = "{}".format(x).split(".")
    if len(y) == 2:
        n = min(n, len(y[-1]))
    return ":.0{:d}f".format(n)


def roundlist(x, max_denominator=1000000):
    x = [
        fractions.Fraction(a).limit_denominator(max_denominator=max_denominator)
        for a in x
    ]
    denoms = [a.denominator for a in x]
    m = functools.reduce(lambda a, b: a * b // fractions.gcd(a, b), denoms)
    return np.array([int(a * m) for a in x])


def weightedsum(values, weights=None):
    if not instance.isarray(values):
        return values
    elif len(values) == 1:
        return values[0]
    elif weights is None or not instance.isarray(weights):
        return np.mean(values)
    else:
        return sum(values * weights) / sum(weights)


def lcm(integers):
    lcm = 1
    for i in integers:
        lcm = lcm * i // math.gcd(lcm, i)
    return lcm

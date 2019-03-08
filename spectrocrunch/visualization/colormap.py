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
import matplotlib.colors as pltcolors

from ..utils import instance


def NormalizedToRGB(x):
    """
    Args:
        x(num|array): data values between 0 and 1
    Retruns:
        r(array):
        g(array):
        b(array):
    """
    x, f = instance.asarrayf(x)
    x = np.round(x*(2**24-1)).astype(int)

    red = f(x & 255)
    green = f((x >> 8) & 255)
    blue = f((x >> 16) & 255)
    return red, green, blue


class LambdaColormap(pltcolors.Colormap):

    def __init__(self, name, func, N=256):
        self._func = func
        super(LambdaColormap, self).__init__(name, N=N)

    def __call__(self, x, alpha=None, bytes=False):
        """
        Args:
            x(num|array): data values between 0 and 
            alpha(Optional(num)): scalar between 0 and 1
            bytes(Optional(bool)): as byte (0-255) or float (0-1)

        Returns:
            tuple: RGBA
        """
        r, g, b = self._func(x)
        if not bytes:
            r = r/255.
            g = g/255.
            b = b/255.
        if instance.isarray(r):
            return np.stack([r, g, b, np.full(r.shape, alpha)], axis=-1)
        else:
            return r, g, b, alpha


def RGBcolormap():
    return LambdaColormap("RGB", NormalizedToRGB, N=2**24-1)


def Linearcolormap(name, a, b, alpha=None):
    a = pltcolors.to_rgba(a, alpha=alpha)
    b = pltcolors.to_rgba(b, alpha=alpha)
    cdict = {'red':  [(0.0, a[0], a[0]),
                      (1.0, b[0], b[0])],
             'green': [(0.0, a[1], a[1]),
                       (1.0, b[1], b[1])],
             'blue': [(0.0, a[2], a[2]),
                      (1.0, b[2], b[2])],
             'alpha': [(0.0, a[3], a[3]),
                       (1.0, b[3], b[3])]}
    return pltcolors.ListedColormap(name, cdict)

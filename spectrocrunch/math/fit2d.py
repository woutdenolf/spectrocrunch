# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
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
import scipy.optimize
import warnings


def gaussian(x, y, x0, y0, sx, sy, rho, A):
    num = (x-x0)**2/sx**2 - 2*rho/(sx*sy)*(x-x0)*(y-y0) + (y-y0)**2/sy**2
    denom = 2*(1-rho**2)
    return A/(2*np.pi*sx*sy*np.sqrt(1-rho**2))*np.exp(-num/denom)


def errorf_gaussian(p, x, y, data):
    x0, y0, sx, sy, rho, A = tuple(p)
    return np.ravel(gaussian(x, y, x0, y0, sx, sy, rho, A)-data)


def guess_gaussian(x, y, data):
    y0i, x0i = np.unravel_index(np.argmax(data), data.shape)
    y0 = y[y0i, 0]
    x0 = x[0, x0i]

    xv = x[y0i, :]-x0
    yv = data[y0i, :]
    sx = np.sqrt(abs(xv**2*yv).sum()/yv.sum())
    xv = y[:, x0i]-y0
    yv = data[:, x0i]
    sy = np.sqrt(abs(xv**2*yv).sum()/yv.sum())
    rho = 0.

    A = data[y0, x0]*2*np.pi*sx*sy*np.sqrt(1-rho**2)

    return np.array([x0, y0, sx, sy, rho, A], dtype=np.float32)


def fitgaussian(x, y, data):
    guess = guess_gaussian(x, y, data)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p, success = scipy.optimize.leastsq(
            errorf_gaussian, guess, args=(x, y, data))
        success = success > 0 and success < 5

    return p, success

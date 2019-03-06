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
from .fit2d import fitgaussian as fitgaussian2d
from .fit1d import fitgaussian as fitgaussian1d


def fmax(data):
    if data.size in data.shape:
        return np.nanargmax(data)
    else:
        return np.array(np.unravel_index(np.nanargmax(data), data.shape))


def fmin(data):
    if data.size in data.shape:
        return np.nanargmin(data)
    else:
        return np.array(np.unravel_index(np.nanargmin(data), data.shape))


def _centroid(data):
    if data.size in data.shape:
        x = np.arange(data.size).reshape(data.shape)
        return (x*data).sum()/data.sum()
    else:
        ny, nx = np.shape(data)
        y, x = np.indices((ny, nx))
        cx = np.sum(x*data)/np.sum(data)
        cy = np.sum(y*data)/np.sum(data)
        return np.array((cx, cy))


def _fit(data):
    if data.size in data.shape:
        x = np.arange(data.size)
        p, success = fitgaussian1d(x, data)
        if success:
            ret = p[0]
    else:
        y, x = np.indices(data.shape)
        p, success = fitgaussian2d(x, y, data)
        if success:
            ret = p[[0, 1]]
    if success:
        return ret
    else:
        return fmax(data)


def foptimize(data, proc, threshold=0.9):
    shift = fmax(data)
    thres = threshold * np.nanmax(data)

    if data.size in data.shape:
        shifta = shift
        shiftb = shift
        while (data.flat[shifta] > thres) and (data.flat[shiftb] > thres):
            shifta -= 1
            shiftb += 1
            if shifta < 0 or shiftb >= data.size:
                shifta += 1
                shiftb -= 1
                break
        if shifta != shiftb:
            shift = proc(data.flat[shifta:shiftb+1])+shifta
    else:
        off = 0
        s = data.shape
        while (data[shift[0]-off:shift[0]+off+1, shift[1]-off:shift[1]+off+1] < thres).sum(dtype=int) == 0:
            off += 1
            if shift[0] < off or shift[1] < off or shift[0]+off >= s[0] or shift[1]+off >= s[1]:
                off -= 1
                break

        if off != 0:
            shift = shift + \
                proc(data[shift[0]-off:shift[0]+off+1,
                          shift[1]-off:shift[1]+off+1])-off

    return shift


def fcentroid(data):
    return foptimize(data, _centroid)


def fgaussmax(data):
    return foptimize(data, _fit)

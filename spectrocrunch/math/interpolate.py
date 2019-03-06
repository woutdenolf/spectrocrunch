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
import scipy.interpolate
from scipy import array
import logging

logger = logging.getLogger(__name__)


def extrap1d(interpolator):
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return array(map(pointwise, array(xs)))

    return ufunclike


def interp1d_floor(xs, ys):

    def pointwise(x):
        ind = np.where(x >= xs)[0]
        if ind.size == 0:
            return np.nan
        else:
            return ys[ind[-1]]

    def ufunclike(x):
        try:
            return array(map(pointwise, array(x)))
        except TypeError:
            return pointwise(x)

    return ufunclike


def interp1d_ceil(xs, ys):

    def pointwise(x):
        ind = np.where(x <= xs)[0]
        if ind.size == 0:
            return np.nan
        else:
            return ys[ind[0]]

    def ufunclike(x):
        try:
            return array(map(pointwise, array(x)))
        except TypeError:
            return pointwise(x)

    return ufunclike


def interpolate_regular(data, axold, axnew, cval=np.nan, degree=1, asgrid=True):
    """Regular data (not necessarily even spaced)

    Args:
        data(array): nD-array
        axold(tuple(array)): 1D-array (regular but not necessary evenly spaced)
        axnew(tuple(array)): 1D-array
        asgrid(Optional(True))

    Returns:
        data(array): (data.size,nD) or mD-array (asgrid)
    """
    ndim = data.ndim
    if len(axold) != ndim or len(axnew) != ndim:
        raise ValueError('Data and axes dimensions must be the same')

    post = None
    args = axnew
    kwargs = {}
    if ndim == 1:
        if degree > 3:
            logger.warning('interpolation degree is capped at 3 (cubic)')
        kind = ['nearest', 'linear', 'quadratic', 'cubic'][min(degree, 3)]
        # nearest==zero, linear==slinear ???
        interp = scipy.interpolate.interp1d(axold[0], data, kind=kind, assume_sorted=False,
                                            fill_value=cval, bounds_error=False)
    elif ndim == 2 and degree > 0:
        interp = scipy.interpolate.RectBivariateSpline(
            axold[0], axold[1], data, kx=degree, ky=degree)
        kwargs['grid'] = asgrid
    else:
        if degree == 0:
            method = 'nearest'
        else:
            if degree > 1:
                logger.warning('interpolation degree is capped at 1 (linear)')
            method = 'linear'
        interp = scipy.interpolate.RegularGridInterpolator(axold, data, method=method,
                                                           fill_value=cval, bounds_error=False)
        if asgrid:
            shape = tuple([len(ax) for ax in axnew])
            def post(x): return x.reshape(shape)
            axnew = np.meshgrid(*axnew, indexing='ij')
            axnew = tuple([ax.flat for ax in axnew])
        args = (np.array(list(zip(*axnew))),)

    ret = interp(*args, **kwargs)
    if post:
        ret = post(ret)
    return ret


def _ravel(ax):
    if isinstance(ax, np.ndarray):
        return ax.ravel()
    else:
        return np.array(ax)


def _ravel_reshape(ax, i, ndim):
    ind = [np.newaxis]*ndim
    ind[i] = slice(None)
    return _ravel(ax)[tuple(ind)]


def interpolate_irregular(data, axold, axnew, cval=np.nan, degree=1, asgrid=True):
    """Irregular data

    Args:
        data(array): nD-array
        axold(tuple(array)): nD-coordinates
        axnew(tuple(array)): 1D-coordinates
        asgrid(Optional(True))

    Returns:
        grid(array): (data.size,nD) or mD-array (asgrid)
    """
    ndim = data.ndim
    if len(axold) != ndim or len(axnew) != ndim:
        raise ValueError('Data and axes dimensions must be the same')

    if ndim == 1:
        kind = ['nearest', 'linear', 'quadratic', 'cubic'][min(degree, 3)]
        # nearest==zero, linear==slinear ???
        interp = scipy.interpolate.interp1d(axold[0], data, kind=kind, assume_sorted=False,
                                            fill_value=cval, bounds_error=False)
        return interp(axnew[0])
    else:
        if ndim == 2:
            method = ['nearest', 'linear', None, 'cubic'][min(degree, 3)]
            if method is None:
                method = 'cubic'
        else:
            if degree == 0:
                method = 'nearest'
            else:
                method = 'linear'
        axold = tuple([_ravel(ax) for ax in axold])
        if asgrid:
            axnew = tuple([_ravel_reshape(ax, i, ndim)
                           for i, ax in enumerate(axnew)])
        else:
            axnew = tuple([_ravel(ax) for ax in axnew])
        return scipy.interpolate.griddata(axold, data.ravel(), axnew, method=method,
                                          fill_value=cval)

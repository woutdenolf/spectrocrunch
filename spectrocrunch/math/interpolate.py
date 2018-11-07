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

def interpolate_ndgrid(data,axold,axnew,cval=np.nan,degree=1,asgrid=True):
    """
    Args:
        data(array): the grid
        axold(list(array)): grid coordinates (regular but not necessary evenly spaced)
        axnew(list(array)): 
    """
    ndim = data.ndim
    if len(axold)!=ndim or len(axnew)!=ndim:
        raise ValueError('Data and axes dimensions must be the same')
            
    post = None
    args = axnew
    kwargs = {}
    if ndim==1:
        kind = ['nearest','linear','quadratic','cubic'][min(degree,3)]
        # nearest==zero, linear==slinear ???
        interp = scipy.interpolate.interp1d(axold[0],data,kind=kind,assume_sorted=False,
                                       fill_value=cval,bounds_error=False)
    elif ndim==2 and degree>0:
        interp = scipy.interpolate.RectBivariateSpline(axold[0],axold[1],data,kx=degree,ky=degree)
        kwargs['grid'] = asgrid
    else:
        if degree==0:
            method = 'nearest'
        else:
            method = 'linear'
        interp = scipy.interpolate.RegularGridInterpolator(axold,data,method=method,
                                                fill_value=cval,bounds_error=False)
        if asgrid:
            shape = tuple([len(ax) for ax in axnew])
            post = lambda x:x.reshape(shape)
            axnew = np.meshgrid(*axnew, indexing='ij')
            axnew = [ax.flat for ax in axnew]
        args = (np.array(zip(*axnew)),)

    ret = interp(*args,**kwargs)
    if post:
        ret = post(ret)
    return ret

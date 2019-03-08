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

from abc import ABCMeta, abstractmethod, abstractproperty
from future.utils import with_metaclass
from contextlib import contextmanager
import numpy as np

from ..math.interpolate import interpolate_regular


class RegularGrid(object):

    DEFAULT_STACKDIM = 0

    def __init__(self, axes, data, stackdim=None):
        self.axes = axes
        self.data = data
        if stackdim is None:
            self.stackdim = self.DEFAULT_STACKDIM
        else:
            self.stackdim = stackdim

    @property
    def shape(self):
        return tuple([ax.size for ax in self.axes])

    def __len__(self):
        return self.axes[0].size

    @property
    def ndim(self):
        return len(self.axes)

    @property
    def size(self):
        return np.prod(self.shape)

    def __repr__(self):
        return str(self.axes)

    def __str__(self):
        return str(self.axes)

    def sliceinfo(self, slicedim):
        """
        Args:
            slicedim(int)
        Returns:
            shape(tuple): shape after slicing
            indexgen(callable): slice index generator
        """
        if slicedim < 0:
            slicedim += self.ndim
        maxdim = self.ndim-1

        if slicedim < 0 or slicedim > maxdim:
            raise ValueError(
                'Slice dimension should be between 0 and {}'.format(slicedim))

        if slicedim == 0:
            def indexgen(i): return (i, Ellipsis)
            shape = self.shape[1:]

        elif slicedim == maxdim:
            def indexgen(i): return (Ellipsis, i)
            shape = self.shape[:-1]
        else:
            a = (slice(None),)*slicedim
            b = (slice(None),)*(maxdim-slicedim)
            def indexgen(i): return a+(i,)+b
            shape = self.shape[:slicedim]+self.shape[slicedim+1:]

        return shape, indexgen

    @contextmanager
    def open(self, **openparams):
        yield self.data

    def locate(self, *ordinates):
        return tuple([ax.locate(x) for x, ax in zip(ordinates, self.axes)])

    def __getitem__(self, index):
        with self.open(mode='r') as data:
            try:
                return data[index]
            except ValueError as e:
                raise IndexError(e)

    def __setitem__(self, index, value):
        with self.open() as data:
            data[index] = value

    def __iadd__(self, value):
        with self.open() as data:
            data[index] += value

    def __isub__(self, value):
        with self.open() as data:
            data[index] -= value

    def __imul__(self, value):
        with self.open() as data:
            data[index] *= value

    def __idiv__(self, value):
        with self.open() as data:
            data[index] /= value

    @property
    def dtype(self):
        with self.open(mode='r') as data:
            return data.dtype

    @property
    def values(self):
        with self.open(mode='r') as data:
            return data

    def interpolate(self, *axes, **kwargs):
        ndim = self.ndim
        if len(axes) != ndim:
            raise ValueError('Expected {} dimensions'.format(ndim))
        axes = [axold.interpolate(axnew)
                for axold, axnew in zip(self.axes, axes)]
        axold, axnew = zip(*axes)
        return interpolate_regular(self, axold, axnew, **kwargs)

    @property
    def stack_dims(self):
        return [self.stackdim]

    @property
    def other_dims(self):
        return [i for i in range(self.ndim) if i not in self.stack_dims]

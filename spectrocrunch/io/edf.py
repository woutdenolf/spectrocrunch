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
import os
import fabio


class edfmemmap():
    """Access edf data with memmaps (cannot handle certain things like compression)
    """

    def __init__(self, filename, mode='r'):
        #f = fabio.edfimage.EdfImage(filename)
        f = fabio.open(filename)

        self.dtype = f.bytecode
        self.shape = (f.dim2, f.dim1)
        self.ndim = len(self.shape)
        offset = f._frames[f.currentframe].start
        self.mdata = np.memmap(filename, dtype=self.dtype,
                               offset=offset, shape=self.shape, order='C')

        if f.swap_needed():
            self.mdata.byteswap(True)

    @property
    def data(self):
        return self.mdata

    def __getitem__(self, index):
        return self.data[index]


class edfimage():
    """Access edf data with fabio
    """

    def __init__(self, filename, mode='r'):
        #self.f = fabio.edfimage.EdfImage(filename)
        self.f = fabio.open(filename)
        self.ndim = 2

    @property
    def dtype(self):
        return self.f.bytecode

    @property
    def shape(self):
        return (self.f.dim2, self.f.dim1)

    @property
    def data(self):
        return self.f.data

    def __getitem__(self, index):
        return self.data[index]

    @property
    def header(self):
        return self.f.header


def saveedf(filename, data, header, overwrite=False):
    exists = os.path.exists(filename)
    if exists:
        if overwrite:
            os.remove(filename)
        else:
            raise IOError("File exists (overwrite=False): {}".format(filename))
    fabio.edfimage.EdfImage(data=data, header=header).write(filename)

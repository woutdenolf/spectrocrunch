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

import os
from contextlib import contextmanager
import numpy as np

from . import axis
from .regulargrid import RegularGrid
from ..io import spec
from ..io.edf import edfimage
from ..instruments import configuration
from ..utils import units

class EDFRegularGrid(RegularGrid):
    
    def __init__(self,filenames,labels=None,**kwargs):
        instrument = configuration.getinstrument(**kwargs)
        parser = spec.edfheader_parser(units=instrument.units,
                        compensationmotors=instrument.compensationmotors,
                        axesnamemap=instrument.imagemotors,
                        **instrument.edfheaderkeys)
        header = edfimage(filenames[0]).header
        info = parser(header)
        axes = info['axes']
        if len(axes) != 2:
            raise RuntimeError('Axes not specified in header')
        
        stackdim = 0
        if labels is None:
            labels = ['.'.join(os.path.basename(f).split('.')[:-1])
                      for f in filenames]
        axnew = axis.factory(labels,type='nominal')
        axes.insert(stackdim,axnew)
        
        self.filenames = filenames
        super(EDFRegularGrid,self).__init__(axes,None,stackdim=stackdim)

    @contextmanager
    def open(self,**openparams):
        yield np.stack([edfimage(f,**openparams).data for f in self.filenames],axis=self.stackdim)
        
    @property
    def dtype(self):
        return edfimage(self.filenames[0]).dtype

    @property
    def values(self):
        with self.open(mode='r') as data:
            return data
    
    def __setitem__(self,index,value):
        raise NotImplementedError()


class XIARegularGrid(RegularGrid):
    
    def __init__(self,xiastack,stacklabel='energylabel',**kwargs):
        instrument = configuration.getinstrument(**kwargs)
        parser = spec.edfheader_parser(units=instrument.units,
                        compensationmotors=instrument.compensationmotors,
                        axesnamemap=instrument.imagemotors,
                        **instrument.edfheaderkeys)
        
        axes = None
        for i,xiaimage in enumerate(xiastack):
            header = xiaimage.header(source=instrument.metadata)
            info = parser(header)
            if axes:
                axes[0][i] = info[stacklabel]
            else:
                s = xiastack.dshape
                stackvalue = info[stacklabel].magnitude
                stackunit = info[stacklabel].units
                detectors = [int(det) for det in xiastack.detectors_used]
                axes = [np.full(s[0],stackvalue),
                        info['axes'][0],info['axes'][1],
                        axis.arange(s[3]),
                        axis.factory(detectors)]
        axes[0] = axis.factory(units.Quantity(axes[0],stackunit))
        
        super(XIARegularGrid,self).__init__(axes,xiastack,stackdim=0)

    def __setitem__(self,index,value):
        raise NotImplementedError()

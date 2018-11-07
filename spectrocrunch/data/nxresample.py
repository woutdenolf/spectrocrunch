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
from scipy import interpolate

from . import nxregulargrid
from . import nxtask
from ..math.interpolate import interpolate_ndgrid
from ..utils import units
from . import axis

class Task(nxregulargrid.Task):
    
    def _parameters_defaults(self):
        super(Task,self)._parameters_defaults()
        self._required_parameters('encoders')
        parameters = self.parameters
        parameters['cval'] = parameters.get('cval',np.nan)
        parameters['kind'] = parameters.get('kind','linear')

    def _prepare_process(self):
        super(Task,self)._prepare_process()
        parameters = self.parameters
        
        encdict = parameters['encoders']
        signals = {sig.name:sig for sig in self.grid.signals}
        encoders = []
        
        axes = self.grid.axes
        axes.pop(self.grid.stackdim)

        for axisdim,ax in enumerate(axes):
            encoder = encdict.get(ax.name,None)
            add = None
            if encoder:
                # ax should have an encoder
                signal = signals.get(encoder['counter'],None)
                if signal:
                    # ax does indeed have an encoder
                    # resolution: steps/unit
                    # offset: steps (encoder position when motor poxition is zero)
                    resolution = units.Quantity(encoder['resolution'],units=1/ax.units)
                    offset = units.Quantity(encoder['offset'],units='dimensionless')
                    encoder = axis.factory(ax.values*resolution+offset)
                    add = self._axis_from_encoder(encoder,axisdim,signal)
            encoders.append(add)

    def _axis_from_encoder(self,axis,axisdim,signal):
        """
        Args:
            axis(Axis): in steps
            axisdim(num)
            signal(Path)
        """
        
        with signal.open() as dset:
            shape = dset.shape
            ndim = dset.ndim
            ind = [slice(None)]*ndim
            
            start_expected = axis.start
            end_expected = axis.end
            ind[axisdim] = 0
            if end_expected>start_expected:
                start_measured = np.min(dset[tuple(ind)])
            else:
                start_measured = np.max(dset[tuple(ind)])
            ind[axisdim] = -1
            if end_expected>start_expected:
                end_measured = np.max(dset[tuple(ind)])
            else:
                end_measured = np.min(dset[tuple(ind)])

            print start_expected,end_expected 
            print start_measured,end_measured
            
            
            
            

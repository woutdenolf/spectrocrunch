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
import logging

from . import nxregulargrid
from . import nxtask
from ..math.interpolate import interpolate_ndgrid
from ..utils import units
from . import axis

logger = logging.getLogger(__name__)

class Task(nxregulargrid.Task):
    
    def _parameters_defaults(self):
        super(Task,self)._parameters_defaults()
        self._required_parameters('encoders')
        parameters = self.parameters
        parameters['cval'] = parameters.get('cval',np.nan)
        parameters['degree'] = parameters.get('degree',1)
        parameters['crop'] = parameters.get('crop',False)
        
    def _parameters_filter(self):
        return super(Task,self)._parameters_filter()+['encoders','cval','degree','crop']

    def _prepare_process(self):
        super(Task,self)._prepare_process()
        parameters = self.parameters
        
        encoders = parameters['encoders']
        signals = {sig.name:sig for sig in self.grid.signals}
        lst = []
        axes = self.grid.axes
        axes.pop(self.grid.stackdim)
        for axisdim,ax in enumerate(axes):
            encoder = encoders.get(ax.name,None)
            add = ax,None
            if encoder:
                # axis should have an encoder
                signal = signals.get(encoder['counter'],None)
                if signal:
                    # axis does indeed have an encoder
                    add = self._encoder(ax,axisdim,signal,encoder['resolution'],
                                        offset=encoder.get('offset',None))
            lst.append(add)

        self.axes,self.encoders = zip(*lst)

    def _process_axes(self):
        return self.axes

    def _encoder(self,position,axisdim,signal,resolution,offset=None):
        """
        Args:
            position(Axis): in motor units
            axisdim(num)
            signal(Path)
            resolution(num): in steps/units
            offset(Optional(num)): in encoder steps
        Returns:
            position(Axis): in motor units
            encoder(Axis): in encoder steps
        """
        
        # Encoder(units) = axis(units)*resolution(steps/units) + offset(units)
        resolution = units.Quantity(resolution,units=1/ax.units)
        if offset is None:
            encoder = axis.factory(position.values*resolution)
        else:
            encoder = axis.factory(position.values*resolution+offset)
        
        crop = self.parameters['crop']
        
        with signal.open() as dset:
            shape = dset.shape
            ndim = dset.ndim
            ind = [slice(None)]*ndim

            # Get encoder start/end
            start_expected = encoder.start
            end_expected = encoder.end
            increasing = end_expected>start_expected
            fmin = lambda ind: np.min(dset[tuple(ind)])
            fmax = lambda ind: np.max(dset[tuple(ind)])
            if crop:
                getlim = fmax,fmin
            else:
                getlim = fmin,fmax

            ind[axisdim] = 0
            if increasing:
                start_measured = getlim[0](ind)
            else:
                start_measured = getlim[1](ind)
            ind[axisdim] = -1
            if increasing:
                end_measured = getlim[1](ind)
            else:
                end_measured = getlim[0](ind)

            # Correct expected encoder position (with offset)
            if offset is None:
                m_expected = (start_expected+end_expected)/2.
                m_measured = (start_measured+end_measured)/2.
                offset = m_measured-m_expected
                encoder = axis.factory(encoder.values+offset)
                start_expected = encoder.start
                end_expected = encoder.end

            # Expand/crop axis to capture the measured positions
            if hasattr(encoder,'stepsize'):
                start_stepsize = encoder.stepsize
                end_stepsize = start_stepsize
            elif hasattr(encoder,'stepsizes'):
                start_stepsize = encoder.stepsizes[0]
                end_stepsize = encoder.stepsizes[1]
            else:
                start_stepsize = encoder.values[1]-encoder.values[0]
                end_stepsize = encoder.values[-2]-encoder.values[-1]
            
            x = encoder.magnitude
            nstepsadd = (start_expected-start_measured)/start_stepsize
            nstepsadd = int(np.ceil(nstepsadd.to('dimensionless').magnitude))
            if nstepsadd>0:
                add = start_expected - np.arange(1,nstepsadd+1)[::-1]*start_stepsize
                x = np.append(add.magnitude,x)
            elif nstepsadd<0:
                x = x[-nstepsadd:]
            
            nstepsadd = (end_measured-end_expected)/start_stepsize
            nstepsadd = int(np.ceil(nstepsadd.to('dimensionless').magnitude))
            if nstepsadd>0:
                add = end_expected + np.arange(1,nstepsadd+1)*end_stepsize
                x = np.append(x,add.magnitude)
            elif nstepsadd<0:
                x = x[:nstepsadd]

            if len(x)>0:
                encoder = axis.factory(x)
                position = self._new_axis((encoder.values-offset)/resolution,ax)
                return position,encoder
            else:
                logger.warning('Skip resampling along axis {}'.format(repr(position.name)))
                return position,None

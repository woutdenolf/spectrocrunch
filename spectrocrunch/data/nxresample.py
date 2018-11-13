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
import logging
import contextlib
import itertools
try:
    from contextlib import ExitStack
except ImportError:
    from contextlib2 import ExitStack
    
from . import nxregulargrid
from . import nxtask
from ..math.interpolate import interpolate_irregular
from ..utils import units
from . import axis

logger = logging.getLogger(__name__)

class Task(nxregulargrid.Task):
    
    def _parameters_defaults(self):
        super(Task,self)._parameters_defaults()
        self._required_parameters('encoders')
        parameters = self.parameters
        parameters['sliced'] = False
        parameters['cval'] = parameters.get('cval',np.nan)
        parameters['degree'] = parameters.get('degree',1)
        parameters['crop'] = parameters.get('crop',False)

    def _parameters_filter(self):
        return super(Task,self)._parameters_filter()+['encoders','cval','degree','crop']

    def _prepare_process(self):
        super(Task,self)._prepare_process()
        parameters = self.parameters
        
        # New axes with associated encoder axis and encoder signal (None if missing)
        encoders = parameters['encoders']
        signals = {sig.name:sig for sig in self.grid.signals}
        lst = []
        for axisdim,ax in enumerate(self.signal_axes):
            encdict = encoders.get(ax.name,None)
            encoder,signal,offset = None,None,0
            if encdict:
                # axis should have an encoder
                signal = signals.get(encdict['counter'],None)
                if signal:
                    # axis does indeed have an encoder
                    ax,encoder,offset = self._encoder(ax,axisdim,signal,encdict['resolution'],
                                        offset=encdict.get('offset',None))
            lst.append((ax,encoder,signal,offset))
        self.axes,self.encoders,self.encoder_signals,self.offsets = zip(*lst)

        # Prepare interpolation
        self.interp_kwargs = {'cval':parameters['cval'],
                              'degree':parameters['degree'],
                              'asgrid':True}
        self.resampled_shape = tuple([s if enc is None else len(enc) for s,enc in zip(self.signalin_shape,self.encoders)])
        
    def _process_axes(self):
        return self.axes

    @contextlib.contextmanager
    def _context_encoders(self):
        with ExitStack() as stack:
            self.encoder_datasets = [stack.enter_context(signal.open())
                                     for signal in self.encoder_signals if signal is not None]
            self.encoder_dest = [enc.magnitude for enc in self.encoders if enc is not None]
            yield
    
    @property
    def signalout_shape(self):
        return self.resampled_shape
    
    def _execute_grid(self):
        """
        Returns:
            nxfs._NXprocess
        """
        # Create new axes (if needed)
        axes = self._create_axes(self._process_axes())

        # Create new signals
        with self._context_encoders():
            indices = list(self._process_indices())
            
            results = self.nxprocess.results
            for signalin in self.grid.signals:
                # Create new NXdata if needed
                nxdata = results[signalin.parent.name]
                bnew = not nxdata.exists
                if bnew:
                    nxdata = results.nxdata(name=signalin.parent.name)

                with signalin.open() as dsetin:
                    # Calculate new signal from old signal
                    signalout = nxdata.add_signal(signalin.name,shape=self.signalout_shape,
                                                  dtype=self.signalout_dtype,chunks=True)
                    with signalout.open() as dsetout:
                        for index in indices:
                            self.indexin = index
                            dsetout[index] = self._process_data(dsetin[index])
                    
                if bnew: 
                    nxdata.set_axes(*axes)
            
            results['encoder_offset'].write(data=self.offsets)
            
    def _process_indices(self):
        shape = self.signalin_shape
        indices = [range(shape[i]) if enc is None else [slice(None)]
                   for i,enc in enumerate(self.encoders)]
        return itertools.product(*indices)

    def _process_data(self,data):
        axold = [dset[self.indexin] for dset in self.encoder_datasets]
        return interpolate_irregular(data,axold,self.encoder_dest,**self.interp_kwargs)

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
            offset(Optional(num)): 
        """
        
        # Encoder(units) = axis(units)*resolution(steps/units) + offset(units)
        resolution = units.Quantity(resolution,units=1/position.units)
        if offset is None:
            encoder = axis.factory(position.values*resolution)
        else:
            encoder = axis.factory(position.values*resolution+offset)
        
        crop = self.parameters['crop']
        
        with signal.open() as dset:
            shape = dset.shape
            ndim = dset.ndim
            ind = [slice(None)]*ndim

            # Encoder: expected start/end
            start_expected = encoder.start
            end_expected = encoder.end
            
            # Encoder: measured start/end (depending on crop)
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

            # Encoder (and corresponding motor positions) at which to interpolate the data
            encoder = encoder.newlimits(start_measured,end_measured)
            if encoder is None:
                logger.warning('Skip resampling along axis {} (dimension was reduced to zero)'.format(repr(position.name)))
                return position,None,offset
            else:
                position = self._new_axis((encoder.values-offset)/resolution,position)
                return position,encoder,offset.to('dimensionless').magnitude

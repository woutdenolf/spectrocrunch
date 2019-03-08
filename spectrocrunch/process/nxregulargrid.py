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
from abc import abstractmethod
import re
import logging

from . import nxprocess
from . import h5regulargrid
from . import axis
from ..utils import instance
from ..utils import units

logger = logging.getLogger(__name__)


class Task(nxprocess.Task):

    DEFAULT_STACKDIM = 0

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        parameters = self.parameters
        parameters['skip'] = parameters.get('skip', [])
        parameters['sliced'] = parameters.get('sliced', False)
        parameters['stackdim'] = parameters.get(
            'stackdim', self.DEFAULT_STACKDIM)

    def _parameters_filter(self):
        return super(Task, self)._parameters_filter()+['sliced', 'stackdim', 'skip']

    def _execute(self):
        """
        Returns:
            nxfs._NXprocess | None
        """
        if len(self.dependencies) != 1:
            raise RuntimeError(
                'nxregulargrid.Task can only depend on exactly one previous task')
        logger.info('Skip signals: {}'.format(self.parameters['skip']))
        self.grid = h5regulargrid.NXRegularGrid(self.previous_outputs[0])
        self._prepare_process()
        self._execute_grid()
        self._sort()

    def _sort(self):
        it = self.temp_nxresults.iter_is_nxclass('NXdata')
        previous_results = self.previous_outputs[0].results
        for nxdata in it:
            if nxdata.islink:
                continue
            nxdataprev = previous_results[nxdata.name]
            if nxdataprev.exists:
                nxdata.sort_signals(other=nxdataprev)

    @property
    def reference_signal_index(self):
        reference = self.parameters.get('reference', None)
        if reference:
            ax = self.grid.axes[self.grid.stackdim]
            ind = np.array([s.path.endswith(reference) for s in ax])
            if ind.any():
                return np.nonzero(ind)[0][-1]
            else:
                raise ValueError(
                    'Reference "{}" not present in {}'.format(reference, ax))
        else:
            return 0

    @property
    def reference_signal(self):
        ax = self.grid.axes[self.grid.stackdim][self.reference_signal_index]
        return h5regulargrid.NXSignalRegularGrid(ax,
            stackdim=self.parameters['stackdim'])

    @property
    def positioners(self):
        return self.temp_nxprocess.positioners()

    def _execute_grid(self):
        """
        Returns:
            nxfs._NXprocess
        """
        # Create new axes (if needed)
        axes = self._create_axes(self._process_axes())

        # Create new signals
        for signalin in self.grid.signals:
            if not self._prepare_signal(signalin):
                continue

            # Create new NXdata if needed
            nxdata = self.temp_nxresults[signalin.parent.name]
            bnew = not nxdata.exists
            if bnew:
                nxdata = self.temp_nxresults.nxdata(name=signalin.parent.name)

            with signalin.open() as dsetin:
                # Calculate new signal from old signal
                if self.sliced:
                    signalout = nxdata.add_signal(signalin.name, shape=self.signalout_shape,
                                                  dtype=self.signalout_dtype, chunks=True)
                    with signalout.open() as dsetout:
                        for i in range(self.signal_nslices):
                            self.indexin[self.signal_stackdim] = i
                            self.indexout[self.signal_stackdim] = i
                            data = self._process_data(
                                dsetin[tuple(self.indexin)])
                            dsetout[tuple(self.indexout)] = data
                else:
                    data = self._process_data(dsetin[tuple(self.indexin)])
                    signalout = nxdata.add_signal(name=signalin.name, data=data,
                                                  chunks=True)

            if bnew:
                nxdata.set_axes(*axes)

    def _prepare_process(self):
        n = self.grid.ndim-1
        self.indexin = [slice(None)]*n
        self.indexout = [slice(None)]*n
        self.skipfuncs = [self._rematch_func(
            redict) for redict in self.parameters['skip']]

    def _skip(self, signal):
        for func in self.skipfuncs:
            if func(signal):
                return True
        return False

    def _prepare_signal(self, signal):
        skip = self._skip(signal)
        if skip:
            logger.info('skip {}'.format(signal.name))
        else:
            logger.info('process {}'.format(signal.name))
        return not skip

    def _process_axes(self):
        return self.signal_axes

    def _create_axes(self, axes):
        """
        Args:
            list(Axis)
        Returns:
            list(str)
        """
        positioners = self.positioners
        for ax in axes:
            if ax.name not in positioners:
                positioners.add_axis(ax.name, ax.values, title=ax.title)
        return [ax.name for ax in axes]

    def _process_data(self, data):
        return data

    @property
    def signal_stackdim(self):
        return self.parameters["stackdim"]

    @property
    def sliced(self):
        return self.parameters["sliced"]

    @property
    def signal_nslices(self):
        return self.signalin_shape[self.signal_stackdim]

    @property
    def dtype_process(self):
        return np.float32(0)

    @property
    def signalout_dtype(self):
        x = np.array(0, dtype=self.grid.dtype)*self.dtype_process
        return x.dtype

    @property
    def signalin_shape(self):
        shape = list(self.grid.shape)
        shape.pop(self.grid.stackdim)
        return tuple(shape)

    @property
    def signalout_shape(self):
        return self.signalin_shape

    @property
    def signal_axes(self):
        axes = list(self.grid.axes)
        axes.pop(self.grid.stackdim)
        return axes

    def _new_axis(self, newvalues, axold):
        #name = '{}_{}'.format(axold.name,self.output.name)
        name = axold.name
        if not isinstance(newvalues, axis.Axis):
            newvalues = units.Quantity(newvalues, units=axold.units)
        return axis.factory(newvalues, name=name, title=axold.title)

    @staticmethod
    def _rematch_func(redict):
        method = redict.get('method', 'regex')
        if method == 'regexparent':
            return lambda signal: re.match(redict['pattern'], signal.parent.name)
        elif method == 'regex':
            return lambda signal: re.match(redict['pattern'], signal.name)
        else:
            return lambda signal: False

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

from . import nxregulargrid
from . import nxtask
from ..utils import instance
from ..io import fs

class Task(nxregulargrid.Task):
    
    def _parameters_defaults(self):
        super(Task,self)._parameters_defaults()
        self._required_parameters('alignmethod','reference')
        parameters = self.parameters
        parameters['refimageindex'] = parameters.get('refimageindex',-1)
        parameters['crop'] = parameters.get('crop',False)
        parameters['roi'] = parameters.get('roi',None)
        parameters['plot'] = parameters.get('plot',False)
        parameters['pad'] = not parameters['crop']
        parameters['onraw'] = True
        
        alignmethod = parameters['alignmethod']
        if alignmethod=="sift":
            from ..align.alignSift import alignSift as alignclass
        elif alignmethod=="elastix":
            from ..align.alignElastix import alignElastix as alignclass
        elif alignmethod=="fft":
            from ..align.alignFFT import alignFFT as alignclass
        elif alignmethod=="min":
            from ..align.alignSimple import alignMin as alignclass
        elif alignmethod=="centroid":
            from ..align.alignSimple import alignCentroid as alignclass
        elif alignmethod=="max":
            from ..align.alignSimple import alignMax as alignclass
        else:
            raise nxtask.ParameterError('Unknown alignmethod {}'.format(repr(alignmethod)))
        self.alignclass = alignclass
    
    def _parameters_filter(self):
        return super(Task,self)._parameters_filter()+['alignmethod','reference','refimageindex','crop','roi','plot','pad','onraw']
        
    def _prepare_reference(self):
        parameters = self.parameters
        refimageindex = parameters['refimageindex']
        refgrid = self.reference_signal
        if instance.isstring(refimageindex):
            if refimageindex=='first':
                refimageindex = 0
            elif refimageindex=='last':
                refimageindex = -1
            elif refimageindex=='middle':
                refimageindex = refgrid.shape[refgrid.stackdim]//2
            else:
                raise nxtask.ParameterError('Unknown refimageindex {}'.format(repr(refimageindex)))
        elif instance.isinteger(refimageindex):
            pass
        else: # fraction
            refimageindex = int(round(refimageindex*refgrid.shape[refgrid.stackdim]))
        parameters['refimageindex'] = refimageindex
        
    def _prepare_process(self):
        # Reference stack and image for alignment
        self._prepare_reference()
        
        # New nxprocess (return when already exists)
        self.process = self._temp_process()

        # Output signal paths
        results = self.process.results
        self.signalsout = []
        for signalin in self.grid.signals:
            nxdata = results[signalin.parent.name]
            if not nxdata.exists:
                nxdata = results.nxdata(name=signalin.parent.name)
            self.signalsout.append(nxdata[signalin.name])
    
    def _process_axes(self,o):
        axes1 = self.signal_axes
        axes2 = o.transform_axes([ax.values for ax in axes1])
        out = []
        for ax1,ax2 in zip(axes1,axes2):
            ax2 = self._new_axis(ax2,ax1)
            if ax1==ax2:
                out.append(ax1)
            else:
                out.append(ax2)
        return out
    
    def _execute_grid(self):
        # Align image stacks
        with self.grid.open_signals() as datasets:
            parameters = self.parameters
            kwargs = {k:parameters[k] for k in ['stackdim','plot']}
            with self.signalsout[0].h5open() as fout:
                signalsout = [sig.path for sig in self.signalsout]
                o = self.alignclass(datasets,None,fout,signalsout,"",**kwargs)
                kwargs = {k:parameters[k] for k in ['onraw','pad','crop','roi','refimageindex']}
                o.align(self.reference_signal_index,**kwargs)
                axes = self._process_axes(o)
                axes = self._create_axes(axes)
                self.process.results['change-of-frame'].write(data=o.absolute_cofs())

        # Set NXdata signal and axes attributes
        for signalout in self.signalsout:
            nxdata = signalout.parent
            with signalout.open() as dsetout:
                signalout = nxdata.add_signal(signalout.name)
            if not nxdata.axes: 
                nxdata.set_axes(*axes)

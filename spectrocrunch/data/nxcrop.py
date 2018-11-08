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
from . import nxregulargrid
from . import nxtask
from . import axis

class Task(nxregulargrid.Task):
    
    def _parameters_defaults(self):
        super(Task,self)._parameters_defaults()
        self._required_parameters('reference')
        if all(p not in self.parameters for p in ['roi','nanval']):
            raise nxtask.MissingParameter('Specify either "nanval" or "roi"')

    def _parameters_filter(self):
        return super(Task,self)._parameters_filter()+['roi','nanval','reference']

    def _prepare_process(self):
        if 'nanval' in self.parameters:
            self.roi = self.calccroproi(self.reference_signal)
        elif 'roi' in self.parameters:
            self.roi = self.convertuserroi(self.reference_signal)

        self.cropped_shape = tuple([b-a for a,b in self.roi])
        self.indexin = [slice(a,b) for a,b in self.roi]
        self.indexout = [slice(None)]*len(self.roi)
    
    @property
    def signal_shape(self):
        return self.cropped_shape
    
    def _process_axes(self):
        axes = []
        for ax,(a,b) in zip(self.signal_axes,self.roi):
            if ax.size!=b-a:
                ax = self._new_axis(ax[a:b],ax)
            axes.append(ax)
        return axes
    
    def calccroproi(self,refgrid):
        """Determine crop ROI to remove slices other than stackdim which contain
        only nanval.
        
        Args:
            refgrid(RegularGrid):

        Returns:
            list(2-tuple): dimensions of refgrid
        """
        
        nanval = self.parameters['nanval']
        
        # Mask (True = valid pixel)
        if self.sliced:
            shape,indexgen = refgrid.sliceinfo(self.signal_stackdim)
            mask = np.ones(shape,dtype=np.bool)
            for i in range(refgrid.shape[self.signal_stackdim]):
                img = refgrid[indexgen(i)]
                if nanval is np.nan:
                    mask &= np.logical_not(np.isnan(img))
                else:
                    mask &= img != nanval
        else:
            if nanval is np.nan:
                mask = np.isnan(refgrid.values).sum(axis=self.signal_stackdim)==0
            else:
                mask = (refgrid.values==nanval).sum(axis=self.signal_stackdim)==0
            shape = mask.shape

        imask = -1
        roi = []
        for igrid in range(refgrid.ndim):
            if igrid==self.signal_stackdim:
                iroi = (0,refgrid.shape[self.signal_stackdim])
            else:
                imask += 1
                sumdims = tuple([i for i in range(refgrid.ndim-1) if i!=imask])
                indvalid = np.argwhere(mask.sum(axis=sumdims))[:,0]
                if indvalid.size==0:
                    return None
                iroi = indvalid[0],indvalid[-1]+1
            roi.append(iroi)

        return roi
    
    def convertuserroi(self,refgrid):
        """Determine crop ROI to remove slices other than stackdim which contain
        only nanval.
        
        Args:
            refgrid(RegularGrid):

        Returns:
            list(2-tuple): dimensions of refgrid
        """
        roi = list(self.parameters['roi'])
        stackdim = self.parameters['stackdim']
        roi.insert(stackdim,(0,refgrid.shape[stackdim]))
        return roi
    

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
from . import regulargrid
from . import nxutils
from . import nxtask

class Task(nxtask.Task):
    
    def _parameters_defaults(self):
        super(Task,self)._parameters_defaults()
        parameters = self.parameters
        parameters['sliced'] = parameters.get('sliced',True)
        parameters['stackdim'] = parameters.get('stackdim',nxutils.DEFAULT_STACKDIM)
        if all(p not in parameters for p in ['roi','nanval']):
            raise nxtask.MissingParameter('Specify either "nanval" or "roi"')
        if 'stackdim' not in parameters:
            raise nxtask.MissingParameter('"stackdim" is missing')

    def _execute(self):
        """
        Returns:
            nxfs._NXprocess | None
        """
        self.grid = regulargrid.NXRegularGrid(self.previous)
        refgrid = self.refgrid
        
        if "nanval" in self.parameters:
            roi = self.calccroproi(refgrid)
        elif "roi" in self.parameters:
            roi = convertuserroi(refgrid)
        else:
            roi = None
            
        if roi:
            return self._crop(roi)
        else:
            return None
            
    @property
    def refgrid(self):
        reference = self.parameters["reference"]
        ax = self.grid.axes[0]
        ind = np.array([s.path.endswith(reference) for s in ax])
        if ind.any():
            ind = np.where(ind)[0][-1]
            return regulargrid.NXSignalRegularGrid(ax[ind])
        else:
            raise ValueError('Reference "{}" not present in {}'.format(ref,ax))
        
    def _crop(self,roi):
        """Apply ROI
        
        Args:
            roi(list(2-tuple)):

        Returns:
            nxfs._NXprocess
        """
        # New nxprocess (return when already exists)
        process,notnew = nxutils.next_process([self.grid.nxgroup],self.parameters)
        if notnew:
            return process
        results = process.results
            
        # Create new axes (if needed)
        positioners = self.grid.nxgroup.positioners()
        axes = []
        gaxes = self.grid.axes
        gaxes.pop(self.grid.stackdim)
        for ax,(a,b) in zip(gaxes,roi):
            if ax.size==b-a:
                axes.append(ax.name)
            else:
                name = '{}_{}'.format(ax.name,self.name)
                positioners.add_axis(name,ax[a:b],title=ax.title)
                axes.append(name)

        # Cropped signals
        stackdim = self.parameters["stackdim"]
        shape = tuple([b-a for a,b in roi])
        indexcrop = [slice(a,b) for a,b in roi]
        indexout = [slice(None)]*len(roi)
        nslice = self.grid.shape[stackdim]
        dtype = self.grid.dtype
        sliced = self.parameters["sliced"]
        for signalin in self.grid.signals:
            nxdata = results[signalin.parent.name]
            bnew = not nxdata.exists
            if bnew:
                nxdata = results.nxdata(name=signalin.parent.name)

            with signalin.open() as dsetin:
                if sliced:
                    signalout = nxdata.add_signal(signalin.name,shape=shape,
                                                  dtype=dtype,chunks=True)
                    with signalout.open() as dsetout:
                        for i in range(nslice):
                            indexcrop[stackdim] = i
                            indexout[stackdim] = i
                            dsetout[tuple(indexout)] = dsetin[tuple(indexcrop)]
                else:
                    signalout = nxdata.add_signal(signalin.name,data=dsetin[tuple(indexcrop)],
                                                  chunks=True)
        
            if bnew: 
                nxdata.set_axes(*axes)
            
        return process
    
    def calccroproi(self,refgrid):
        """Determine crop ROI to remove slices other than stackdim which contain
        only nanval.
        
        Args:
            refgrid(RegularGrid):

        Returns:
            tuple or None
        """
        
        stackdim = self.parameters['stackdim']
        sliced = self.parameters['sliced']
        nanval = self.parameters['nanval']
        
        # Mask (True = valid pixel)
        if sliced:
            shape,indexgen = refgrid.sliceinfo(stackdim)
            mask = np.ones(shape,dtype=np.bool)
            for i in range(refgrid.shape[stackdim]):
                img = refgrid[indexgen(i)]
                if nanval is np.nan:
                    mask &= np.logical_not(np.isnan(img))
                else:
                    mask &= img != nanval
        else:
            if nanval is np.nan:
                mask = np.isnan(refgrid.values).sum(axis=stackdim)==0
            else:
                mask = (refgrid.values==nanval).sum(axis=stackdim)==0
            shape = mask.shape

        imask = -1
        roi = []
        for igrid in range(refgrid.ndim):
            if igrid==stackdim:
                iroi = (0,refgrid.shape[stackdim])
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
            list(2-tuple)
        """
        roi = list(self.parameters["roi"])
        stackdim = self.parameters['stackdim']
        roi.insert(stackdim,(0,refgrid.shape[stackdim]))
        return roi
    

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

from ..io import spec
from ..io import nexus
from ..io import edf
from ..patch.pint import ureg
from ..utils import instance
from ..utils import units
from ..utils import listtools
from ..instruments import configuration
from ..process.edfregulargrid import EDFRegularGrid
from ..process.nxresult import regulargriddata

import numpy as np
import warnings
import os
import collections


class Base(object):

    def __init__(self,axis0name=None,axis1name=None,
                 transfo0=None,transfo1=None,**kwargs):
        self.instrument = configuration.getinstrument(kwargs)
        if axis0name is None:
            axis0name = self.instrument.imageaxes[0].upper()
        if axis1name is None:
            axis1name = self.instrument.imageaxes[1].upper()
        self.axis0name = axis0name
        self.axis1name = axis1name
        if transfo0 is None:
            transfo0 = lambda x:x
        if transfo1 is None:
            transfo1 = lambda x:x
        self.transfo0 = transfo0
        self.transfo1 = transfo1

    def motorquantity(self,values,motorname):
        return units.Quantity(values,self.instrument.units[motorname])


class PointBase(Base):

    def setcoordinates(self):
        coord0 = []
        coord1 = []
        lbls = []
        for positions,labels in self:
            p0 = sum(self.motorquantity(x,motname) for motname,x in
                        zip(self.instrument.imagemotors,positions)
                        if self.instrument.imageaxes[0] in motname)
            p1 = sum(self.motorquantity(x,motname) for motname,x in
                        zip(self.instrument.imagemotors,positions)
                        if self.instrument.imageaxes[1] in motname)

            # Append positions and labels
            coord0 += instance.asarray(p0).tolist()
            coord1 += instance.asarray(p1).tolist()
            lbls += labels
        if not coord0:
            raise RuntimeError("No Base found")
        self._coordinates0 = units.asqarray(coord0)
        self._coordinates1 = units.asqarray(coord1)
        self.labels = lbls

    @property
    def coordinates0(self):
        return self.transfo0(self._coordinates0)

    @property
    def coordinates1(self):
        return self.transfo1(self._coordinates1)
        
        
class ImageBase(Base):

    def __init__(self,grid,**kwargs):
        self.instrument = configuration.getinstrument(kwargs)
        self.grid = grid
        axis0name = self.grid.axes[1].name
        axis1name = self.grid.axes[2].name
        if self.transpose:
            axis0name,axis1name = axis1name,axis0name
        kwargs['axis0name'] = kwargs.get('kwargs',axis0name).upper()
        kwargs['axis1name'] = kwargs.get('kwargs',axis1name).upper()
        kwargs['instrument'] = self.instrument
        super(ImageBase,self).__init__(**kwargs)
    
    @property
    def transpose(self):
        ax0 = self.grid.axes[0][0]
        if ax0 not in self.instrument.imageaxes:
            return False
        else:
            return self.grid.axes[0][0] != self.instrument.imageaxes[0]
    
    @property
    def axis0values(self):
        return self.transfo0(self.grid.axes[1].values)

    @property
    def axis1values(self):
        return self.transfo1(self.grid.axes[2].values)

    def displaydata(self,index=None):
        """
        Args:
            index(Optional(list)): e.g. [5,None,1]
            
        Returns:
            data(array): e.g. nrow x ncol x 2
            channels: e.g. [0,None,1]
            labels: e.g. ["group5","group1"]
        """
        if index is None:
            index = range(len(self))
        else:
            index = instance.asarray(index).tolist()
        nimages = len(index)-index.count(None)
        
        data = np.zeros(self.grid.shape[1:]+(nimages,),dtype=self.grid.dtype)
        labels = [""]*nimages
        channels = list(index)
        iout = -1
        for j,ind in enumerate(index):
            if ind is None:
                continue
            iout += 1
            labels[iout] = self.grid.axes[0][ind]
            data[...,iout] = self.grid[ind,...]
            channels[j] = iout
            
        if self.transpose:
            data = np.swapaxes(data,0,1)
        return data,channels,labels

    def interpolate(self,p0,p1):
        ind = list(range(len(p0)))
        data = self.grid.interpolate(None,p0,p1,asgrid=True,degree=1)
        data = data[:,ind,ind]
        result = collections.OrderedDict()
        for label,values in zip(self.grid.axes[0].values,data):
            result[label] = values
        return result
        

class EDFStack(ImageBase):

    def __init__(self,filenames,**kwargs):
        """
        Args:
            filename(str|list(str)): list of edf file names
        """
        if instance.isstring(filenames):
            filenames = [filenames]
        grid = EDFRegularGrid(filenames,instrument=kwargs.get('instrument',None))
        super(EDFStack,self).__init__(grid,**kwargs)
        
        
class NexusStack(ImageBase):

    def __init__(self,nxgroup,**kwargs):
        grid = regulargriddata(nxgroup)
        super(NexusStack,self).__init__(grid,**kwargs)
        
        
class XanesSpec(PointBase):

    def __init__(self,filenames,specnumbers,labels=None,**kwargs):
        """
        Args:
            filenames(str|list(str)): list of spec file names
            specnumbers(list|list(list)): empty list of numbers => all xanes spectra
            labels(Optional(list|list(list))): uses the spec numbers by default
        """
        super(XanesSpec,self).__init__(**kwargs)
        
        if instance.isstring(filenames):
            filenames = [filenames]
        self.filenames = filenames
        self.specnumbers = self.listoflists(specnumbers)
        self.speclabels = self.listoflists(labels)
        self.setcoordinates()

    @staticmethod
    def listoflists(lst):
        if lst is None:
            return [[]]
        if not instance.isarray(lst):
            lst = [lst]
        if lst:
            if not instance.isarray(lst[0]):
                lst = [lst]
        else:
            lst = [lst]
        return lst

    def __iter__(self):
        for filename,numbers,labels in zip(self.filenames,self.specnumbers,self.speclabels):
            # Get motor positions for each number
            f = spec.spec(filename)
            if not numbers:
                numbers = f.extractxanesginfo(keepsum=True,sumingroups=False,keepindividual=False)
                numbers = [k[0] for k in numbers if len(k)==1]
                lbls = []
            if not numbers:
                continue
            positions = zip(*[f.getmotorvalues(nr,self.instrument.imagemotors) for nr in numbers])
            if not labels:
                labels = numbers
            yield positions,labels

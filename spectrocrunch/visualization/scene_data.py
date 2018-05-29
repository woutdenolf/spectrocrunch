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
from .. import ureg
from ..common import instance
from ..common import units
from ..common import listtools
from ..instruments import configuration
from ..h5stacks.get_hdf5_imagestacks import get_hdf5_imagestacks

import numpy as np
import warnings
import os
import scipy.interpolate
import collections

class Coordinates(object):

    def __init__(self,instrument=None,axis0name=None,axis1name=None,transfo0=None,transfo1=None):
        self.instrument = configuration.factory(instrument)
        
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

    def motorpositions(self,values,motorname):
        return units.Quantity(values,self.instrument.units[motorname])


class PointCoordinates(Coordinates):

    def setcoordinates(self):
        coord0 = []
        coord1 = []
        lbls = []
        for positions,labels in self:
            p0 = sum(self.motorpositions(x,motname) for motname,x in zip(self.instrument.imagemotors,positions) if self.instrument.imageaxes[0] in motname)
            p1 = sum(self.motorpositions(x,motname) for motname,x in zip(self.instrument.imagemotors,positions) if self.instrument.imageaxes[1] in motname)

            # Append positions and labels
            coord0.extend(instance.asarray(p0).tolist())
            coord1.extend(instance.asarray(p1).tolist())
            lbls.extend(labels)

        if not coord0:
            raise RuntimeError("No coordinates found")

        self._coordinates0 = units.asqarray(coord0)
        self._coordinates1 = units.asqarray(coord1)
        
        self.labels = lbls

    @property
    def coordinates0(self):
        return self.transfo0(self._coordinates0)

    @property
    def coordinates1(self):
        return self.transfo1(self._coordinates1)
        
        
class ImageCoordinates(Coordinates):

    def __init__(self,labels=None,**kwargs):
        super(ImageCoordinates,self).__init__(**kwargs)
        if labels is None:
            labels = []
        self.labels = labels
    
    def iter_interpolate(self):
        return iter(self)
    
    def setcoordinates(self,axis0values,axis1values,motdim0,motdim1,offsets):
        # Scan positions
        axis0values = self.motorpositions(axis0values,motdim0)
        axis1values = self.motorpositions(axis1values,motdim1)
        
        # Scan offsets
        for mot in self.instrument.compensationmotors[motdim0]:
            off,func = instance.asarrayf(offsets[mot][0])
            if not instance.isnumber(off[0]):
                off = off.astype(float)
            off = func(off)
            axis0values = axis0values+self.motorpositions(off,mot)
            
        for mot in self.instrument.compensationmotors[motdim1]:
            off,func = instance.asarrayf(offsets[mot][0])
            if not instance.isnumber(off[0]):
                off = off.astype(float)
            off = func(off)
            axis1values = axis1values+self.motorpositions(off,mot)
        
        # Make sure that the vertical axis is dim0
        transpose = self.instrument.imageaxes[0] in motdim1
        if transpose:
            axis0values,axis1values = axis1values,axis0values

        self._axis0values = axis0values
        self._axis1values = axis1values
        self.transpose = transpose
        
        self.axis0f = None
        self.axis1f = None
    
    @property
    def axis0values(self):
        return self.transfo0(self._axis0values)

    @property
    def axis1values(self):
        return self.transfo1(self._axis1values)

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

        data = None
        labels = [""]*nimages
        channels = list(index)
        
        iout = -1
        for j,i in enumerate(index):
            if i is None:
                continue
            else:
                iout += 1
            channels[j] = iout
            
            image,label = self[i]
            
            try:
                labels[iout] = self.labels[i]
            except IndexError:
                labels[iout] = label

            if data is None:
                if nimages==1:
                    data = image
                    break
                else:
                    data = np.zeros(image.shape+(nimages,),dtype=image.dtype)
            data[...,iout] = image
     
        if self.transpose:
            data = np.swapaxes(data,0,1)
            
        return data,channels,labels

    def interpolate(self,p0,p1):
        unit0 = self.axis0values.units
        unit1 = self.axis1values.units

        if self.axis0f is None:
            self.axis0f = scipy.interpolate.interp1d(self.axis0values.magnitude,np.arange(len(self.axis0values)),\
                                                    kind='linear',bounds_error=False,fill_value=np.nan)
                                                    
        if self.axis1f is None:
            self.axis1f = scipy.interpolate.interp1d(self.axis1values.magnitude,np.arange(len(self.axis1values)),\
                                                    kind='linear',bounds_error=False,fill_value=np.nan)
                                                    
        p0 = instance.asarray(units.magnitude(p0,unit0))
        p1 = instance.asarray(units.magnitude(p1,unit1))
        result = collections.OrderedDict()

        x0 = self.axis0f(p0)
        x1 = self.axis1f(p1)
        indvalid = np.isfinite(x0) & np.isfinite(x1)
        x0v = np.round(x0[indvalid]).astype(np.int)
        x1v = np.round(x1[indvalid]).astype(np.int)

        values = None
        for data,label in self.iter_interpolate():
            values = np.full(x0.size,np.nan)
            values[indvalid] = data[x0v,x1v]
            result[label] = values
        
        return result
        

class ZapRoiMap(ImageCoordinates):

    def __init__(self,filenames,**kwargs):
        """
        Args:
            filename(str|list(str)): list of edf file names
        """
        super(ZapRoiMap,self).__init__(**kwargs)
        
        if instance.isstring(filenames):
            filenames = [filenames]
        self.filenames = filenames
        self.setmotorinfo()

    def __len__(self):
        return len(self.filenames)
    
    def __getitem__(self,index):
        filename = self.filenames[index]
        label = os.path.splitext(os.path.basename(filename))[0]
        f = edf.edfimage(filename)
        return f.data,label
            
    def setmotorinfo(self):
        # Parse edf header
        o = spec.cmd_parser()
        f = edf.edfimage(self.filenames[0])
        header = f.header
        cmd = header[self.instrument.edfheaderkeys["speccommand"]]
        result = o.parsezapimage(cmd)

        # Extract axes
        if result["name"]!='zapimage':
            raise RuntimeError("Cannot extract zapimage information from \"{}\"".format(cmd))
        axis0values = spec.zapline_values(result["startfast"],result["endfast"],result["npixelsfast"])
        axis1values = spec.ascan_values(result["startslow"],result["endslow"],result["nstepsslow"])

        # Store coordinates
        self.setcoordinates(axis0values,axis1values,result["motfast"],result["motslow"],header)
        
        # Image saved: slow x fast
        self.transpose = not self.transpose


class Nexus(ImageCoordinates):

    def __init__(self,filename,groups,defaultstackindex=None,**kwargs):
        """
        Args:
            filename(str): h5 file name
            groups(list(dict)): {"path":str,"ind":int}
            defaultstackindex(Optional(int)): 
        """
        super(Nexus,self).__init__(**kwargs)
        self.filename = filename
        if not instance.isarray(groups):
            groups = [groups]
        self.groups = groups
        if defaultstackindex is None:
            defaultstackindex = -1
        self.defaultstackindex = defaultstackindex
        self.setmotorinfo()
    
    def __len__(self):
        return len(self.groups)

    def __getitem__(self,index):
        with nexus.File(self.filename,mode="r") as oh5:
            group = self.groups[index]
            label = group["path"].split("/")[-1]
            if group.get("ind",None) is None:
                group["ind"] = self.defaultstackindex
            try:
                data,axes,axesnames = nexus.parse_NXdata(oh5[group["path"]])
            except KeyError:
                raise KeyError("{} not in {}".format(group["path"],filename))
            stackindex = self.stackindex(axesnames)
            if stackindex==0:
                data = data[group["ind"],...]
            elif stackindex==1:
                data = data[:,group["ind"],:]
            else:
                data = data[...,group["ind"]]
    
        return data,label
        
    def iter_interpolate(self):
        stacks,axes,procinfo = get_hdf5_imagestacks(self.filename,self.instrument.h5stackgroups)
        stackindex = self.stackindex([k["name"] for k in axes])
        
        with nexus.File(self.filename,mode="r") as oh5:
            for k1 in stacks:
                for k2 in stacks[k1]:
                    data,axes,axesnames = nexus.parse_NXdata(oh5[stacks[k1][k2]])
                    n = data.shape[stackindex]
                    for i in range(n):
                        if n>1:
                            label = "{} ({})".format(k2,axes[stackindex][i])
                        else:
                            label = k2
                        if stackindex==0:
                            yield data[i,...],label
                        elif stackindex==1:
                            yield data[:,i,:],label
                        else:
                            yield data[...,i],label
                
    def stackindex(self,axesnames):
        stackindex = [name for name in axesnames if name not in self.instrument.units]
        if len(stackindex)!=1:
            raise RuntimeError("Could not determine the stack index: {}".format(axesnames))
        return axesnames.index(stackindex[0])

    def setmotorinfo(self):
        with nexus.File(self.filename,mode="r") as oh5:
            # Get motor info from hdf5
            try:
                ocoord = oh5["stackinfo"]
            except KeyError:
                warnings.warn("\"coordinates\" is deprecated and should be replaced by \"stackinfo\"", DeprecationWarning) 
                ocoord = oh5["coordinates"]
            _,axes,axesnames = nexus.parse_NXdata(oh5[self.groups[0]["path"]])
            
            # Get scanning axes
            b = [name in self.instrument.units for name in axesnames]
            if sum(b)!=2:
                raise RuntimeError("Expected exactly two scanning dimensions: {}".format(axesnames))
            axes = listtools.listadvanced_bool(axes,b)
            axesnames = listtools.listadvanced_bool(axesnames,b)

            # Store coordinates
            self.setcoordinates(axes[0][:],axes[1][:],axesnames[0],axesnames[1],ocoord)
        
        
class XanesSpec(PointCoordinates):

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
        

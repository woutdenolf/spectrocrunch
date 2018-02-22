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

import numpy as np
import warnings
import os

INSTRUMENTINFO = {}
INSTRUMENTINFO["id21"] = {"IMAGEMOTORS" : ["samy","sampy","samz","sampz"],\
                            "IMAGEAXES" : ["z","y"],\
                            "UNITS" : {"samy":ureg.millimeter,\
                                    "samz":ureg.millimeter,\
                                    "samx":ureg.millimeter,\
                                    "sampy":ureg.micrometer,\
                                    "sampz":ureg.micrometer},\
                            "MOTOR_OFFSETS" : {"samy":["sampy"],\
                                      "samz":["sampz"],\
                                      "sampy":["samy"],\
                                      "sampz":["samz"]},\
                            "SPECCMDLABEL":"Title"}


class Coordinates(object):

    def __init__(self,instrument=None):
        self.instrument = instrument

    def __getattr__(self,attr):
        return INSTRUMENTINFO[self.instrument][attr]

    def motorpositions(self,values,motorname):
        return units.Quantity(values,self.UNITS[motorname])


class PointCoordinates(Coordinates):

    def setcoordinates(self):
        coord0 = []
        coord1 = []
        lbls = []
        for positions,labels in self:
            p0 = sum(self.motorpositions(x,motname) for motname,x in zip(self.IMAGEMOTORS,positions) if self.IMAGEAXES[0] in motname)
            p1 = sum(self.motorpositions(x,motname) for motname,x in zip(self.IMAGEMOTORS,positions) if self.IMAGEAXES[1] in motname)

            # Append positions and labels
            coord0.append(p0)
            coord1.append(p1)
            lbls.extend(labels)

        if not coord0:
            raise RuntimeError("No coordinates found")
        
        # There is not extend for pint.Quantity
        self.coordinates0 = units.flatten(coord0)
        self.coordinates1 = units.flatten(coord1)
        self.labels = lbls


class ImageCoordinates(Coordinates):

    def __init__(self,axis0name=None,axis1name=None,labels=None,**kwargs):
        super(ImageCoordinates,self).__init__(**kwargs)

        if axis0name is None:
            axis0name = self.IMAGEAXES[0].upper()
        if axis1name is None:
            axis1name = self.IMAGEAXES[1].upper()
        self.axis0name = axis0name
        self.axis1name = axis1name
        
        if labels is None:
            labels = []
        self.labels = labels
        
    def setcoordinates(self,axis0values,axis1values,motdim0,motdim1,offsets):
        # Scan positions
        axis0values = self.motorpositions(axis0values,motdim0)
        axis1values = self.motorpositions(axis1values,motdim1)
        
        # Scan offsets
        for mot in self.MOTOR_OFFSETS[motdim0]:
            off,func = instance.asarrayf(offsets[mot])
            if not instance.isnumber(off[0]):
                off = off.astype(float)
            off = func(off)
            axis0values = axis0values+self.motorpositions(off,mot)
        for mot in self.MOTOR_OFFSETS[motdim1]:
            off,func = instance.asarrayf(offsets[mot])
            if not instance.isnumber(off[0]):
                off = off.astype(float)
            off = func(off)
            axis1values = axis1values+self.motorpositions(off,mot)
        
        # Make sure that the vertical axis is dim0
        transpose = self.IMAGEAXES[0] in motdim1
        if transpose:
            axis0values,axis1values = axis1values,axis0values

        self.axis0values = axis0values
        self.axis1values = axis1values
        self.transpose = transpose

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
        cmd = header[self.SPECCMDLABEL]
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
        with nexus.File(self.filename) as oh5:
            group = self.groups[index]
            label = group["path"].split("/")[-1]
            if group.get("ind",None) is None:
                group["ind"] = self.defaultstackindex
            
            try:
                data,axes,axesnames = nexus.parse_NXdata(oh5[group["path"]])
            except KeyError:
                raise KeyError("{} not in {}".format(group["path"],filename))
        
            stackindex = [name for name in axesnames if name not in self.UNITS]
            if len(stackindex)!=1:
                raise RuntimeError("Could not determine the stack index: {}".format(axesnames))
            stackindex = stackindex[0]

            if stackindex==0:
                data = data[group["ind"],...]
            elif stackindex==1:
                data = data[:,group["ind"],:]
            else:
                data = data[...,group["ind"]]
    
        return data,label
    
    def setmotorinfo(self):
        with nexus.File(self.filename) as oh5:
            # Get motor info from hdf5
            try:
                ocoord = oh5["stackinfo"]
            except KeyError:
                warnings.warn("\"coordinates\" is deprecated and should be replaced by \"stackinfo\"", DeprecationWarning) 
                ocoord = oh5["coordinates"]
            _,axes,axesnames = nexus.parse_NXdata(oh5[self.groups[0]["path"]])
            
            # Get scanning axes
            b = [name in self.UNITS for name in axesnames]
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
            filename(str|list(str)): list of edf file names
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

            positions = zip(*[f.getmotorvalues(nr,self.IMAGEMOTORS) for nr in numbers])

            if not labels:
                labels = numbers

            yield positions,labels
        

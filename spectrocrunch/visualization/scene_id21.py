# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
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

from . import scene
from ..io import spec
from ..io import nexus
from ..io import edf
from .. import ureg
from ..common import instance
from ..common import units

import numpy as np
import warnings

class MotorCoordinates(object):

    IMAGEMOTORS = ["samy","sampy","samz","sampz"]
    IMAGEAXES = ["z","y"] # dim0, dim1
    
    UNITS = {"samy":ureg.millimeter,\
            "samz":ureg.millimeter,\
            "samx":ureg.millimeter,\
            "sampy":ureg.micrometer,\
            "sampz":ureg.micrometer}
            
    MOTOR_OFFSETS = {"samy":["sampy"],\
              "samz":["sampz"],\
              "sampy":["samy"],\
              "sampz":["samz"]}

    @classmethod
    def defaultinit(cls,kwargs):
        kwargs["axisname0"] = kwargs.get("axisname0",cls.IMAGEAXES[0].upper())
        kwargs["axisname1"] = kwargs.get("axisname1",cls.IMAGEAXES[1].upper())

    @classmethod
    def motorpositions(cls,values,motorname):
        return units.Quantity(values,cls.UNITS[motorname])

    @classmethod
    def coordinates(cls,axis0values,axis1values,motdim0,motdim1,offsets):
        # Scan positions
        axis0values = cls.motorpositions(axis0values,motdim0)
        axis1values = cls.motorpositions(axis1values,motdim1)
        
        # Scan offsets
        for mot in cls.MOTOR_OFFSETS[motdim0]:
            off,func = instance.asarrayf(offsets[mot])
            if not instance.isnumber(off[0]):
                off = off.astype(float)
            off = func(off)
            axis0values = axis0values+cls.motorpositions(off,mot)
        for mot in cls.MOTOR_OFFSETS[motdim1]:
            off,func = instance.asarrayf(offsets[mot])
            if not instance.isnumber(off[0]):
                off = off.astype(float)
            off = func(off)
            axis1values = axis1values+cls.motorpositions(off,mot)
        
        # Make sure that the vertical axis is dim0
        transpose = cls.IMAGEAXES[0] in motdim1
        if transpose:
            axis0values,axis1values = axis1values,axis0values

        return axis0values,axis1values,transpose


class ZapRoiMap(scene.Image,MotorCoordinates):

    def __init__(self,filenames,**kwargs):
        """
        Args:
            filename(str|list(str)): list of edf file names
        """
        self.defaultinit(kwargs)
        img,axis0values,axis1values,labels = self.parseedf(filenames)
        super(ZapRoiMap,self).__init__(img,lim0=axis0values[[0,-1]],lim1=axis1values[[0,-1]],labels=labels,**kwargs)
    
    @classmethod
    def parseedf(cls,filenames):
        if instance.isstring(filenames):
            filenames = [filenames]
        
        for filename in filenames:
            if filename is None:
                continue
            f = edf.edfimage(filename)
            break
        else:
            raise RuntimeError("At least one filename must be given")
             
        axis0values,axis1values,transpose = cls.limits(f.header,"Title")

        if len(filenames)>1:
            tmp = f.data
            img = np.zeros(tmp.shape+(len(filenames),),dtype=tmp.dtype)

            for i,filename in enumerate(filenames):
                if filename is None:
                    continue
                f = edf.edfimage(filename)
                img[...,i] = f.data
            if transpose:
                img = np.swapaxes(img,0,1)
        else:
            img = f.data
            if transpose:
                img = img.T
        
        labels = os.path.splitext(os.path.basename(filenames))[0]
        
        return img,axis0values,axis1values,labels
    
    @classmethod
    def limits(cls,header,cmdlabel):
        o = spec.cmd_parser()
        cmd = header[cmdlabel]
        result = o.parsezapimage(cmd)

        if result["name"]!='zapimage':
            raise RuntimeError("Cannot extract zapimage information from \"{}\"".format(cmd))
        axis0values = spec.zapline_values(result["startfast"],result["endfast"],result["npixelsfast"])
        axis1values = spec.ascan_values(result["startslow"],result["endslow"],result["nstepsslow"])
        
        axis0values,axis1values,transpose = cls.coordinates(axis0values,axis1values,result["motfast"],result["motslow"],header)
        
        # Image saved: slow x fast
        return axis0values,axis1values,not transpose


class Nexus(scene.Image,MotorCoordinates):

    def __init__(self,filename,groups,**kwargs):
        """
        Args:
            filename(str): h5 file name
            groups(dict): stackname:stackindex
        """
        self.defaultinit(kwargs)
        img,axis0values,axis1values,labels = self.parseh5(filename,groups)
        super(Nexus,self).__init__(img,lim0=axis0values,lim1=axis1values,labels=labels,**kwargs)

    @classmethod
    def parseh5(cls,filename,groups):
        if not instance.isarray(groups):
            groups = [groups]
            
        oh5 = nexus.File(filename)
        
        try:
            ocoord = oh5["stackinfo"]
        except KeyError:
            warnings.warn("\"coordinates\" is deprecated and should be replaced by \"stackinfo\"", DeprecationWarning) 
            ocoord = oh5["coordinates"]

        img = None
        labels = []
        for igroup,group in enumerate(groups):
            if group is None:
                labels.append("")
                continue
            else:
                labels.append(group["path"].split("/")[-1])
            try:
                data,axes,axesnames = nexus.parse_NXdata(oh5[group["path"]])
            except KeyError:
                raise KeyError("{} not in {}".format(group["path"],filename))
                
            if img is None:
                stackindex = np.where([a not in cls.UNITS for a in axesnames])[0][0]
                if stackindex==0:
                    data = data[group["ind"],...]
                elif stackindex==1:
                    data = data[:,group["ind"],:]
                else:
                    data = data[...,group["ind"]]
  
                axes.pop(stackindex)
                axesnames.pop(stackindex)

                axis0values,axis1values,transpose = cls.coordinates(axes[0][:],axes[1][:],axesnames[0],axesnames[1],ocoord)
                if transpose:
                    data = data.T
                
                if len(groups)==1:
                    img = data
                else:
                    img = np.zeros(data.shape+(len(groups),),dtype=data.dtype)
                    img[...,igroup] = data
            else:
                if stackindex==0:
                    data = data[group["ind"],...]
                elif stackindex==1:
                    data = data[:,group["ind"],:]
                else:
                    data = data[...,group["ind"]]
                if transpose:
                    data = data.T
                    
                img[...,igroup] = data

        oh5.close()
        
        return img,axis0values,axis1values,labels


class XanesSpec(scene.Text,MotorCoordinates):

    def __init__(self,filenames,specnumbers,labels=None,**kwargs):
        """
        Args:
            filename(str|list(str)): list of edf file names
            specnumbers(list|list(list)): empty list of numbers => all xanes spectra
        """
        self.defaultinit(kwargs)
        coord0,coord1,labels = self.parsespec(filenames,specnumbers,labels=labels)
        super(XanesSpec,self).__init__(coord0,coord1,labels=labels,**kwargs)

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

    @classmethod
    def parsespec(cls,filenames,specnumbers,labels=None):
        
        if instance.isstring(filenames):
            filenames = [filenames]
        
        specnumbers = cls.listoflists(specnumbers)
        speclabels = cls.listoflists(labels)

        coord0 = []
        coord1 = []
        labels = []
        for filename,numbers,lbls in zip(filenames,specnumbers,speclabels):
            # Get motor positions for each number
            f = spec.spec(filename)
            if not numbers:
                numbers = f.extractxanesginfo(keepsum=True,sumingroups=False,keepindividual=False)
                numbers = [k[0] for k in numbers if len(k)==1]
                lbls = []
            if not numbers:
                continue
            
            # Parse motor positions
            positions = zip(*[f.getmotorvalues(nr,cls.IMAGEMOTORS) for nr in numbers])
            p0 = sum(cls.motorpositions(x,motname) for motname,x in zip(cls.IMAGEMOTORS,positions) if cls.IMAGEAXES[0] in motname)
            p1 = sum(cls.motorpositions(x,motname) for motname,x in zip(cls.IMAGEMOTORS,positions) if cls.IMAGEAXES[1] in motname)

            # Append positions and labels
            coord0.append(p0)
            coord1.append(p1)
            if lbls:
                labels.extend(lbls)
            else:
                labels.extend(numbers)
        
        if not coord0:
            raise RuntimeError("No XANES scans found in spec files: {}".format(filenames))
        
        # There is not extend for pint.Quantity
        coord0 = units.flatten(coord0)
        coord1 = units.flatten(coord1)

        return coord0,coord1,labels
        

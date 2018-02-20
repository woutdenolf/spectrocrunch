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

import numpy as np
import warnings

class MotorCoordinates(object):

    UNITS = {"samy":ureg.millimeter,\
            "samz":ureg.millimeter,\
            "samx":ureg.millimeter,\
            "sampy":ureg.micrometer,\
            "sampz":ureg.micrometer}
    MOTOFF = {"samy":["sampy"],\
              "samz":["sampz"],\
              "sampy":["samy"],\
              "sampz":["samz"]}
              
    VERTAXIS = "z"
    
    @classmethod
    def coordinates(cls,dim0values,dim1values,motdim0,motdim1,offsets):
        # Scan positions
        dim0values = ureg.Quantity(dim0values,cls.UNITS[motdim0])
        dim1values = ureg.Quantity(dim1values,cls.UNITS[motdim1])
        
        # Scan offsets
        for mot in cls.MOTOFF[motdim0]:
            off,func = instance.asarrayf(offsets[mot])
            if not instance.isnumber(off[0]):
                off = off.astype(float)
            off = func(off)
            dim0values = dim0values+ureg.Quantity(off,cls.UNITS[mot])
        for mot in cls.MOTOFF[motdim1]:
            off,func = instance.asarrayf(offsets[mot])
            if not instance.isnumber(off[0]):
                off = off.astype(float)
            off = func(off)
            dim1values = dim1values+ureg.Quantity(off,cls.UNITS[mot])
        
        # Make sure that the vertical axis is dim0
        transpose = cls.VERTAXIS in motdim1
        if transpose:
            dim0values,dim1values = dim1values,dim0values

        return dim0values,dim1values,transpose


class ZapRoiMap(scene.Image,MotorCoordinates):

    def __init__(self,filenames,**kwargs):
        """
        Args:
            filename(str|list(str)): list of edf file names
        """
        
        if instance.isstring(filenames):
            filenames = [filenames]
        
        for filename in filenames:
            if filename is None:
                continue
            f = edf.edfimage(filename)
            break
        else:
            raise RuntimeError("At least one filename must be given")
             
        dim0values,dim1values,transpose = self.limits(f.header,"Title")

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
            
        super(ZapRoiMap,self).__init__(img,lim0=dim0values[[0,-1]],lim1=dim1values[[0,-1]],dim0name="Z",dim1name="Y",**kwargs)
    
    @classmethod
    def limits(cls,header,cmdlabel):
        o = spec.cmd_parser()
        cmd = header[cmdlabel]
        result = o.parsezapimage(cmd)

        if result["name"]!='zapimage':
            raise RuntimeError("Cannot extract zapimage information from \"{}\"".format(cmd))
        dim0values = spec.zapline_values(result["startfast"],result["endfast"],result["npixelsfast"])
        dim1values = spec.ascan_values(result["startslow"],result["endslow"],result["nstepsslow"])
        
        dim0values,dim1values,transpose = cls.coordinates(dim0values,dim1values,result["motfast"],result["motslow"],header)
        
        # Image saved: slow x fast
        return dim0values,dim1values,not transpose


class Nexus(scene.Image,MotorCoordinates):

    def __init__(self,filename,groups,**kwargs):
        """
        Args:
            filename(str): h5 file name
            groups(dict): stackname:stackindex
        """
        if not instance.isarray(groups):
            groups = [groups]
            
        oh5 = nexus.File(filename)
        
        try:
            ocoord = oh5["stackinfo"]
        except KeyError:
            warnings.warn("\"coordinates\" is deprecated and should be replaced by \"stackinfo\"", DeprecationWarning) 
            ocoord = oh5["coordinates"]

        img = None
        for igroup,group in enumerate(groups):
            if group is None:
                continue
            try:
                data,axes,axesnames = nexus.parse_NXdata(oh5[group["path"]])
            except KeyError:
                raise KeyError("{} not in {}".format(group["path"],filename))
                
            if img is None:
            
                stackindex = np.where([a not in self.UNITS for a in axesnames])[0][0]
                if stackindex==0:
                    data = data[group["ind"],...]
                elif stackindex==1:
                    data = data[:,group["ind"],:]
                else:
                    data = data[...,group["ind"]]
  
                axes.pop(stackindex)
                axesnames.pop(stackindex)

                dim0values,dim1values,transpose = self.coordinates(axes[0][:],axes[1][:],axesnames[0],axesnames[1],ocoord)
                if transpose:
                    data = data.T
                
                if len(groups)==1:
                    img = data
                else:
                    img = np.zeros(data.shape+(len(groups),),dtype=data.dtype)
                    img[...,0] = data
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

        super(Nexus,self).__init__(img,lim0=dim0values,lim1=dim1values,dim0name="Z",dim1name="Y",**kwargs)
        

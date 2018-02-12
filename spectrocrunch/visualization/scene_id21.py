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
from .. import ureg
from ..common import instance

import fabio
import numpy as np
import warnings

class ZapRoiMap(scene.Image):

    UNITS = {"samy":ureg.millimeter,"samz":ureg.millimeter,"sampy":ureg.micrometer,"sampz":ureg.micrometer}
    MOTOFF = {"samy":["sampy"],\
              "samz":["sampz"],\
              "sampy":["samy"],\
              "sampz":["samz"]}

    def __init__(self,filenames):
        
        if instance.isstring(filenames):
            filenames = [filenames]
        
        for filename in filenames:
            if filename is None:
                continue
            f = fabio.open(filename)
            break

        header = f.header
        lim0,lim1,transpose = self.limits(header,"Title")

        if len(filenames)>1:
            tmp = f.data
            img = np.zeros(tmp.shape+(3,),dtype=tmp.dtype)
            for i,filename in enumerate(filenames):
                if filename is None:
                    continue
                f = fabio.open(filename)
                img[...,i] = f.data
            if transpose:
                img = np.swapaxes(img,0,1)
        else:
            img = f.data
            if transpose:
                img = img.T
            
        super(ZapRoiMap,self).__init__(img,lim0=lim0,lim1=lim1,dim0name="Z",dim1name="Y")
    
    @classmethod
    def limits(self,header,cmdlabel):
        o = spec.cmd_parser()
        cmd = header[cmdlabel]
        result = o.parsezapimage(cmd)

        if result["name"]!='zapimage':
            raise RuntimeError("Cannot extract zapimage information from \"{}\"".format(cmd))
        fastvalues = spec.zapline_values(result["startfast"],result["endfast"],result["npixelsfast"])
        slowvalues = spec.ascan_values(result["startslow"],result["endslow"],result["nstepsslow"])
        
        fastvalues = ureg.Quantity(fastvalues,self.UNITS[result["motfast"]])
        slowvalues = ureg.Quantity(slowvalues,self.UNITS[result["motslow"]])
        
        for mot in self.MOTOFF[result["motfast"]]:
            fastvalues = fastvalues+ureg.Quantity(float(header[mot]),self.UNITS[mot])
        for mot in self.MOTOFF[result["motslow"]]:
            slowvalues = slowvalues+ureg.Quantity(float(header[mot]),self.UNITS[mot])
        
        # Image saved: slow x fast
        transpose = "z" in result["motfast"]
        if transpose:
            lim0 = fastvalues
            lim1 = slowvalues
        else:
            lim0 = slowvalues
            lim1 = fastvalues

        print lim0[[0,-1]]
        print lim1[[0,-1]]
        
        return lim0,lim1,transpose

class Nexus(scene.Image):

    UNITS = {"samy":ureg.millimeter,\
            "samz":ureg.millimeter,\
            "samx":ureg.millimeter,\
            "sampy":ureg.micrometer,\
            "sampz":ureg.micrometer}
    MOTOFF = {"samy":["sampy"],\
              "samz":["sampz"],\
              "sampy":["samy"],\
              "sampz":["samz"]}

    def __init__(self,filename,groups,originzero=False):

        if not instance.isarray(groups):
            groups = [groups]
            
        oh5 = nexus.File(filename)
        
        try:
            ocoord = oh5["stackinfo"]
        except KeyError:
            warnings.warn("\"coordinates\" is deprecated and should be replaced by \"stackinfo\"", DeprecationWarning) 
            ocoord = oh5["coordinates"]
        ocoord = {a:ureg.Quantity(instance.asarray(ocoord[a].value),self.UNITS[a]) for a in ocoord if a in self.UNITS}

        img = None

        for igroup,group in enumerate(groups):
            if group is None:
                continue
            data,axes,axesnames = nexus.parse_NXdata(oh5[group["path"]])

            if img is None:
            
                stackindex = np.where([a not in self.UNITS for a in axesnames])[0][0]
                if stackindex==0:
                    data = data[group["ind"],...]
                elif stackindex==1:
                    data = data[:,group["ind"],:]
                else:
                    data = data[...,group["ind"]]
  
                del axes[stackindex]
                del axesnames[stackindex]

                dim0values = ureg.Quantity(axes[0][:],self.UNITS[axesnames[0]])
                dim1values = ureg.Quantity(axes[1][:],self.UNITS[axesnames[1]])
                dim0mot,dim1mot = axesnames
                
                for mot in self.MOTOFF[dim0mot]:
                    add = ocoord[mot]
                    if len(add)==1:
                        add = add[0]
                    else:
                        add = add[group["ind"]]
                    dim0values = dim0values+add
                for mot in self.MOTOFF[dim1mot]:
                    add = ocoord[mot]
                    if len(add)==1:
                        add = add[0]
                    else:
                        add = add[group["ind"]]
                    dim1values = dim1values+add

                dim1values = dim1values[[0,-1]]
                dim0values = dim0values[[0,-1]]
                if originzero:
                    dim1values -= np.min(dim1values)
                    dim0values -= np.min(dim0values)
                
                transpose = "y" in dim0mot
                if transpose:
                    lim0 = dim1values
                    lim1 = dim0values
                else:
                    lim0 = dim0values
                    lim1 = dim1values
            
                if transpose:
                    data = data.T

                
                if len(groups)==1:
                    img = data
                else:
                    img = np.zeros(data.shape+(3,),dtype=data.dtype)
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

        super(Nexus,self).__init__(img,lim0=lim0,lim1=lim1,dim0name="Z",dim1name="Y")
        

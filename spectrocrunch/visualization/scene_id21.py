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
from .. import ureg
from ..common import instance

import fabio
import numpy as np

class ZapRoiMap(scene.Image):

    UNITS = {"samy":ureg.millimeter,"samz":ureg.millimeter,"sampy":ureg.micrometer,"sampz":ureg.micrometer}
    MOTOFF = {"samy":["sampy"],\
              "samz":["sampz"],\
              "sampy":["samy"],\
              "sampz":["samz"]}

    def __init__(self,filenames):
        
        if instance.isstring(filenames):
            filenames = [filenames]
        
        f = fabio.open(filenames[0])
        
        header = f.header
        lim0,lim1,transpose = self.limits(header,"Title")
        
        if len(filenames)>1:
            tmp = f.data
            img = np.zeros(tmp.shape+(3,),dtype=tmp.dtype)
            img[...,0] = tmp
            for i in range(1,min(len(filenames),3)):
                f = fabio.open(filenames[i])
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
        
        transpose = "z" in result["motfast"]
        if transpose:
            lim0 = fastvalues
            lim1 = slowvalues
        else:
            lim0 = slowvalues
            lim1 = fastvalues
        
        return lim0,lim1,transpose
        

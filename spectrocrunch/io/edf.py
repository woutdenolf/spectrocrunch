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

from . import spec

import numpy as np

import fabio

def get2Dscancoordinates(header,scanlabel="",fastlabel="",slowlabel=""):
    """Get scan coordinates from header

    Args:
        header(dict): edf header
        scanlabel(optional(str)): header key which value is the scan command
        fastlabel(optional(str)): header key which value is the fast motor name
        slowlabel(optional(str)): header key which value is the slow motor name

    Returns:
        (dict,dict): {"name":str,"data":array}
    """

    ret = ({"name":"","data":[]},{"name":"","data":[]})

    scaninfo = {"name":""}

    if len(scanlabel)>0:
        if scanlabel in header:
            cmd = header[scanlabel]
            o = spec.cmd_parser()
            scaninfo = o.parse(cmd)

    if scaninfo["name"] != "zapimage":
        if len(fastlabel)>0 and len(slowlabel)>0:
            label = fastlabel
            if label+"_mot" in header:
                scaninfo["motfast"] = str(header[label+"_mot"])
                scaninfo["startfast"] = np.float(header[label+"_start"])
                scaninfo["endfast"] = np.float(header[label+"_end"])
                scaninfo["npixelsfast"] = np.float(header[label+"_nbp"])
            else:
                return ret

            label = slowlabel
            if label+"_mot" in header:
                scaninfo["motslow"] = str(header[label+"_mot"])
                scaninfo["startslow"] = np.float(header[label+"_start"])
                scaninfo["endslow"] = np.float(header[label+"_end"])
                scaninfo["nstepsslow"] = np.float(header[label+"_nbp"])
            else:
                return ret

            scaninfo["name"] = "zapimage"

    if scaninfo["name"] == "zapimage":
        sfast = {"name":scaninfo["motfast"],"data":spec.zapline_values(scaninfo["startfast"],scaninfo["endfast"],scaninfo["npixelsfast"])}
        sslow = {"name":scaninfo["motslow"],"data":spec.ascan_values(scaninfo["startslow"],scaninfo["endslow"],scaninfo["nstepsslow"])}

        ret = (sfast,sslow)

    return ret

class edfmemmap():
    """Access edf data with memmaps (cannot handle certain things like compression)
    """

    def __init__(self,filename,mode='r'):
        #f = fabio.edfimage.EdfImage(filename)
        f = fabio.open(filename)

        self.dtype = f.bytecode
        self.shape = (f.dim2,f.dim1)
        self.ndim = len(self.shape)
        offset = f._frames[f.currentframe].start
        self.mdata = np.memmap(filename,dtype=self.dtype,offset=offset,shape=self.shape,order='C')

        if f.swap_needed():
            self.mdata.byteswap(True)

    @property
    def data(self):
        return self.mdata
        
    def __getitem__(self,index):
        return self.data[index]

class edffabio():
    """Access edf data with fabio
    """

    def __init__(self,filename,mode='r'):
        #self.f = fabio.edfimage.EdfImage(filename)
        self.f = fabio.open(filename)

        self.ndim = 2

    @property
    def dtype(self):
        return self.f.bytecode

    @property
    def shape(self):
        return (self.f.dim2,self.f.dim1)

    @property
    def data(self):
        return self.f.data
        
    def __getitem__(self,index):
        return self.data[index]



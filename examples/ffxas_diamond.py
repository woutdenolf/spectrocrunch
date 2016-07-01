# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
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

import os, sys
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spectrocrunch.align.alignElastix import alignElastix as align
#from spectrocrunch.align.alignSimple import alignMax as align
#from spectrocrunch.align.alignSift import alignSift as align

import h5py

if __name__ == '__main__':
    path = "/data/visitor/ls2497/id21/XANES from Diamond/"
    source = os.path.join(path,"Xanes1Co.h5")
    dest = os.path.join(path,"Xanes1Co.aligned.h5")
    stackdim = 2
    sourcelist = ["/exchange/data"]
    destlist = ["/exchange/data"]
    roi = ((0,-1),(5,50))

    a = 15
    b = 58
    c = 20
    hin = h5py.File(source, "r")
    source = [hin["/exchange/data"][...,a:b]]
    sourcelist = "None"

    os.remove(dest)
    o = align(source,sourcelist,dest,destlist,"",stackdim=stackdim,overwrite=True,plot=True)
    o.align(0, refimageindex=c, onraw = True, pad = False, crop = True, roi = roi)
    
    hout = h5py.File(dest, "a")
    hout["/exchange/energy"] = hin["/exchange/energy"][a:b]
    hin.close()
    hout.close()

    

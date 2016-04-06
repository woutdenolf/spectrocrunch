# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
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

from spectrocrunch.align.alignElastix import alignElastix
from spectrocrunch.align.alignFFT import alignFFT
from spectrocrunch.align.alignSift import alignSift
from spectrocrunch.align.tests.teststack import teststack

import os
import glob
import numpy as np

def showstack(stackData):
    import PyQt4.Qt as qt
    from PyMca5.PyMca import StackBrowser

    app = qt.QApplication([])

    w = StackBrowser.StackBrowser()
    w.setStackDataObject(stackData, index=0)
    w.show()
    
    app.exec_()

def alignexample(t):

    # Source data (several image stacks)
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),"testdata")

    if t=="xrfxanes":
        source = os.path.join(path,t,"h5","rape5_XANESfull.h5")
        stackdim = 2
        sourcelist = ["DSAg","DSCu","DSIt","DSZn"]
        nstack = len(sourcelist)
        refdatasetindex = 1
        refimageindex = 0
        alignclass = alignElastix
    elif t =="fluotomo":
        source = ""
        dataset1 = sorted(glob.glob(os.path.join(path,t,"*.tif")))
        sourcelist = [dataset1]
        nstack = 1
        refdatasetindex = 0
        refimageindex = 0
        alignclass = alignElastix
    elif t == "testdata":
        source,offsets,stackdim = teststack()
        sourcelist = None
        nstack = len(source)
        refdatasetindex = 0
        refimageindex = 0
        alignclass = alignSift
        print(offsets)
    else:
        return

    # Destination
    outputstack = [np.zeros(1,dtype=np.float32)]*nstack

    # Align
    o = alignclass(source,sourcelist,outputstack,None,None,stackdim=stackdim,overwrite=True)

    #o.align_pairwise(refdatasetindex,onraw = False,extend = False)
    #o.align_pairwise(refdatasetindex,onraw = True,extend = False)
    #o.align_pairwise(refdatasetindex,onraw = False,extend = True)
    o.align_pairwise(refdatasetindex,onraw = True,extend = True)
    #o.align_fixed(refdatasetindex,refimageindex)

    if t == "testdata":
        offsets2 = o.offsets
        offsets2[:,0] -= offsets2[0,0]
        offsets2[:,1] -= offsets2[0,1]
        print(offsets2)

    # Show result
    showstack(np.rollaxis(outputstack[0],stackdim,0))

if __name__ == '__main__':

    alignexample("testdata")



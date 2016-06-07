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

import unittest

from ..alignElastix import alignElastix
from ..alignSift import alignSift
from ..alignFFT import alignFFT
from ..alignSimple import alignMin
from ..alignSimple import alignMax
from ..alignSimple import alignCentroid

from .teststack import teststack

import numpy as np

class test_align(unittest.TestCase):

    def getrelativeoffset(self,offsets,refimageindex):
        off = offsets.copy()
        off[:,0] -= off[refimageindex,0]
        off[:,1] -= off[refimageindex,1]
        return off

    def getoffsetrelativetomedian(self,offsets):
        off = offsets.copy()
        off[:,0] -= np.median(off[:,0])
        off[:,1] -= np.median(off[:,1])
        return off

    def test_align(self,alignclass):
        # Prepare dataIO
        inputstack,offsets,stackdim = teststack()
        outputstack = [np.zeros(1,dtype=np.float32)]*len(inputstack)

        # References
        refdatasetindex = 0
        refimageindex = len(inputstack)//2
        offsets_median = self.getoffsetrelativetomedian(offsets)
        offsets_ref = self.getrelativeoffset(offsets,refimageindex)
        offsets_0 = self.getrelativeoffset(offsets,0)
        
        # Prepare alignment
        o = alignclass(inputstack,None,outputstack,None,None,stackdim=stackdim,overwrite=True)

        # Check alignment
        for i in range(4):
            pad = (i & 1)==1
            crop = (i & 2)==2
            # Pairwise: align on aligned
            o.align(refdatasetindex,onraw = False,pad = pad,crop = crop,roi=((1,-1),(1,-1)))
            offsets_0b = self.getrelativeoffset(o.offsets,0)
            np.testing.assert_almost_equal(offsets_0,offsets_0b,decimal=1)
            # Pairwise: align on raw
            o.align(refdatasetindex,onraw = True,pad = pad,crop = crop,roi=((1,-1),(1,-1)))
            offsets_0b = self.getrelativeoffset(o.offsets,0)
            np.testing.assert_almost_equal(offsets_0,offsets_0b,decimal=1)
            # Fixed reference
            o.align(refdatasetindex,refimageindex=refimageindex,pad = pad,crop = crop,roi=((1,-1),(1,-1)))
            offsets_0b = self.getrelativeoffset(o.offsets,0)
            np.testing.assert_almost_equal(offsets_0,offsets_0b,decimal=1)

    def test_elastix(self):
        self.test_align(alignElastix)

    def test_sift(self):
        self.test_align(alignSift)

    def test_fft(self):
        self.test_align(alignFFT)

    def test_min(self):
        self.test_align(alignMin)

    def test_max(self):
        self.test_align(alignMax)

    def test_centroid(self):
        self.test_align(alignCentroid)

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    #testSuite.addTest(test_align("test_min"))
    testSuite.addTest(test_align("test_max"))
    #testSuite.addTest(test_align("test_centroid"))
    testSuite.addTest(test_align("test_fft"))
    testSuite.addTest(test_align("test_sift"))
    testSuite.addTest(test_align("test_elastix"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

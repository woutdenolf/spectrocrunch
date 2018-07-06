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

from ..types import transformationType
from . import helper_teststack

import numpy as np
import matplotlib.pyplot as plt

class test_teststack(unittest.TestCase):

    def test_data(self):
        types = [transformationType.translation, transformationType.rigid, transformationType.similarity, transformationType.affine, transformationType.homography]
        for t in types:
            if t==transformationType.translation:
                lst = [True,False]
            else:
                lst = [False]
    
            for vector in lst:
                for transposed in lst:
                    listofstacks,COFrelative,stackdim = helper_teststack.data(t,vector=vector,transposed=transposed)
                    self.assertIsInstance(listofstacks,list)
                    self.assertIsInstance(listofstacks[0],np.ndarray)
                    self.assertEqual(len(listofstacks[0].shape),3)
                    self.assertTrue(all(s.shape==listofstacks[0].shape for s in listofstacks))

                    if False:
                        for s in listofstacks:
                            for i in range(s.shape[stackdim]):
                                if stackdim==0:
                                    self.plot(s[i,...])
                                elif stackdim==1:
                                    self.plot(s[:,i,:])
                                else:
                                    self.plot(s[...,i])
                            break # show only one

    def plot(self,img):
        plt.figure(1)
        plt.subplot(111)
        plt.imshow(img,origin='lower',interpolation='nearest')
        plt.pause(0.1)

def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_teststack("test_data"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

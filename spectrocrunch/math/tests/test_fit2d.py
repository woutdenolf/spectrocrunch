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

from .. import fit2d

import numpy as np
import pylab

class test_fit2d(unittest.TestCase):

    def test_leastsq(self):
        nx = 501
        ny = 401
        y, x = np.indices((ny,nx))
        x0 = 10
        y0 = ny//2
        sx = nx//4
        sy = ny//4
        rho = 0.5
        A = 1000.
        p1 = np.array([x0,y0,sx,sy,rho,A],dtype=np.float32)
        x0,y0,sx,sy,rho,A = tuple(p1)

        data = fit2d.gaussian(x,y,x0,y0,sx,sy,rho,A)
        #self.plot(data)

        p2,_ = fit2d.fitgaussian(x,y,data)
        np.testing.assert_allclose(p1,p2)

    def plot(self,img):
        pylab.figure(1)
        pylab.subplot(111)
        pylab.imshow(img,origin='lower',interpolation='nearest')
        pylab.pause(0.1)
        raw_input("Press enter to continue...")

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_fit2d("test_leastsq"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

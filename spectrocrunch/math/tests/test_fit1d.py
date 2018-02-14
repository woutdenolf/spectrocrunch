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

import unittest

from .. import fit1d
from ...simulation import noisepropagation

import numpy as np


class test_fit1d(unittest.TestCase):

    def test_lstsq_std(self):
        nx,npa = 10,3
        A = np.random.rand(nx,npa)*10
        
        # None of these work:
        
        # VAR(b) -LSTSQ-> VAR(x) -LINPROP-> VAR(b)
        x = np.random.rand(npa)
        b = np.dot(A,x)
        varb = b # poisson
        varx = fit1d.lstsq_std(A,None,None,vare=varb)
        x = noisepropagation.randomvariable(x,np.sqrt(varx))
        b2 = np.dot(A,x)
        varb2 = noisepropagation.VAR(b2)
        #np.testing.assert_allclose(varb,varb2)
        
        # VAR(x) -LINPROP-> VAR(b) -LSTSQ-> VAR(x)
        x = noisepropagation.poisson(np.random.rand(npa))
        b2 = np.dot(A,x)
        varb2 = noisepropagation.VAR(b2)
        varx2 = fit1d.lstsq_std(A,None,None,vare=varb2)
        #np.testing.assert_allclose(noisepropagation.VAR(x),varx2)
        
        # x -LINPROP-> VAR(b)
        # VAR(x) -LINPROP-> VAR(b)
        x = np.random.rand(npa)
        varx = np.random.rand(npa)**2
        x = noisepropagation.randomvariable(x,varx)
        varb = noisepropagation.VAR(np.dot(A,x))
        varb2 = np.dot(A*A,varx)
        #np.testing.assert_allclose(varb,varb2)
        
    def test_leastsq(self):
        nx = 501
        x = np.arange(nx)
        x0 = 10
        sx = nx//4
        A = 1000.
        p1 = np.array([x0,sx,A],dtype=np.float32)
        x0,sx,A = tuple(p1)

        data = fit1d.gaussian(x,x0,sx,A)
        #import matplotlib.pyplot as plt
        #plt.plot(data)
        #plt.show()
        
        p2,_ = fit1d.fitgaussian(x,data)
        np.testing.assert_allclose(p1,p2)

    def test_linfit(self):
        x = np.asarray([1.47,1.50,1.52,1.55,1.57,1.60,1.63,1.65,1.68,1.70,1.73,1.75,1.78,1.80,1.83])
        y = np.asarray([52.21,53.12,54.48,55.84,57.20,58.57,59.93,61.29,63.11,64.47,66.28,68.10,69.92,72.19,74.46])
        m = 61.272
        b = -39.062
        vm = 3.1539
        vb = 8.63185
        (m1,b1),(em1,eb1) = fit1d.linfit(x,y,errors=True)
        (m2,b2),(em2,eb2) = fit1d.linfit2(x,y,errors=True)
        vm1 = em1*em1
        vm2 = em2*em2
        vb1 = eb1*eb1
        vb2 = eb2*eb2
        np.testing.assert_array_almost_equal([m,b],[m1,b1],decimal=3)
        np.testing.assert_array_almost_equal([m,b],[m2,b2],decimal=3)
        np.testing.assert_array_almost_equal([vm,vm],[vm1,vm2],decimal=3)
        np.testing.assert_array_almost_equal([vb,vb],[vb1,vb2],decimal=4)
        np.testing.assert_allclose([m1,b1],[m2,b2])
        
        m1,em1 = fit1d.linfit_zerointercept(x,y,errors=True)
        m2,em2 = fit1d.linfit_zerointercept2(x,y,errors=True)
        np.testing.assert_allclose(m1,m2)
        np.testing.assert_almost_equal(em1,em2,decimal=3) # em2 maybe not correct???

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_fit1d("test_leastsq"))
    testSuite.addTest(test_fit1d("test_linfit"))
    testSuite.addTest(test_fit1d("test_lstsq_std"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

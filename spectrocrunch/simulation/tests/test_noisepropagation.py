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

from .. import noisepropagation

import numpy as np

import itertools

class test_noisepropagation(unittest.TestCase):

    def test_repeat(self):
        X = noisepropagation.randomvariable([100,200],[100**0.5,200**0.5])
        
        N = 10
        Y1 = noisepropagation.repeat(N,X)
        Y2 = noisepropagation.compound(noisepropagation.randomvariable(N,0),X)
        Y2 = np.squeeze(Y2)
        
        np.testing.assert_allclose(noisepropagation.E(Y1),noisepropagation.E(Y2))
        np.testing.assert_allclose(noisepropagation.S(Y1),noisepropagation.S(Y2))
        np.testing.assert_allclose(noisepropagation.SNR(X),noisepropagation.SNR(N*X))
        np.testing.assert_allclose(noisepropagation.SNR(Y1),noisepropagation.SNR(X)*N**0.5)
        
    def test_compound(self):
    
        N0 = noisepropagation.randomvariable([100.,200],[100**0.5,200**0.5])
        EN0 = noisepropagation.E(N0)
        
        g1 = noisepropagation.randomvariable(10.,1.3)
        g2 = noisepropagation.randomvariable(3.,0.2)
        g3 = noisepropagation.randomvariable(7.,1.5)
        g4 = noisepropagation.randomvariable(6.,1.7)
        
        Eg1 = noisepropagation.E(g1)
        Eg2 = noisepropagation.E(g2)
        Eg3 = noisepropagation.E(g3)
        Eg4 = noisepropagation.E(g4)
                
        N1 = noisepropagation.compound(N0,g1)[0,:]
        N2 = noisepropagation.compound(N1,g2)[0,:]
        N3 = noisepropagation.compound(N2,g3)[0,:]
        N4 = noisepropagation.compound(N3,g4)[0,:]
        
        a = Eg1*Eg2*Eg3*Eg4
        
        b = noisepropagation.VAR(g1) * (Eg2*Eg3*Eg4)**2 +\
            noisepropagation.VAR(g2) * (Eg1)*(Eg3*Eg4)**2 +\
            noisepropagation.VAR(g3) * (Eg1*Eg2)*(Eg4)**2 +\
            noisepropagation.VAR(g4) * (Eg1*Eg2*Eg3) \
            
        np.testing.assert_allclose(noisepropagation.E(N4) , a*EN0)
        np.testing.assert_allclose(noisepropagation.VAR(N4), b*EN0 + a**2 * noisepropagation.VAR(N0))
        
        np.testing.assert_allclose(noisepropagation.RVAR(N4), noisepropagation.RVAR(N0) +\
                                                             (noisepropagation.RVAR(g1) +\
                                                              noisepropagation.RVAR(g2)/Eg1 +\
                                                              noisepropagation.RVAR(g3)/(Eg1*Eg2) +\
                                                              noisepropagation.RVAR(g4)/(Eg1*Eg2*Eg3)
                                                             )/EN0  )
                                                             
    def test_bernouilli(self):
        # Bernouilli processes can be multiplied
        g1 = noisepropagation.bernouilli(0.5)
        g2 = noisepropagation.bernouilli(0.4)
        g3 = noisepropagation.bernouilli(noisepropagation.E(g1)*noisepropagation.E(g2))
        
        N = noisepropagation.randomvariable(100.,12.3)
        
        N1 = noisepropagation.compound(noisepropagation.compound(N,g1),g2)
        N2 = noisepropagation.compound(N,g3)
        self.assertEqual(noisepropagation.E(N1),noisepropagation.E(N2))
        self.assertEqual(noisepropagation.S(N1),noisepropagation.S(N2))
        
        # Bernouilli compounding of Poisson gives Poisson
        N = noisepropagation.poisson(100)
        N1 = noisepropagation.compound(N,g1)
        N2 = noisepropagation.poisson(noisepropagation.E(N)*noisepropagation.E(g1))
        self.assertEqual(noisepropagation.E(N1),noisepropagation.E(N2))
        self.assertEqual(noisepropagation.S(N1),noisepropagation.S(N2))
        
    def test_poisson(self):
        # Poisson processes cannot be multiplied
        g1 = noisepropagation.poisson(10)
        g2 = noisepropagation.poisson(15)
        g3 = noisepropagation.poisson(noisepropagation.E(g1)*noisepropagation.E(g2))
        
        N = noisepropagation.randomvariable(100.,12.3)
        
        N1 = noisepropagation.compound(noisepropagation.compound(N,g1),g2)
        N2 = noisepropagation.compound(N,g3)
        self.assertEqual(noisepropagation.E(N1),noisepropagation.E(N2))
        self.assertNotEqual(noisepropagation.S(N1),noisepropagation.S(N2))
     
        # Poisson compounding of Poisson does not gives Poisson
        N = noisepropagation.poisson(100)
        N1 = noisepropagation.compound(N,g1)
        N2 = noisepropagation.poisson(noisepropagation.E(N)*noisepropagation.E(g1))
        self.assertEqual(noisepropagation.E(N1),noisepropagation.E(N2))
        self.assertNotEqual(noisepropagation.S(N1),noisepropagation.S(N2))
        
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_noisepropagation("test_repeat"))
    testSuite.addTest(test_noisepropagation("test_compound"))
    testSuite.addTest(test_noisepropagation("test_bernouilli"))
    testSuite.addTest(test_noisepropagation("test_poisson"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

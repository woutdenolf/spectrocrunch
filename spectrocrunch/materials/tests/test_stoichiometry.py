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

from ..import stoichiometry

import numpy as np

class test_stoichiometry(unittest.TestCase):

    def test_fraction(self):
        n = np.asarray([0.1,0.6,0.3])
        MM = np.asarray([5.6,3.5,7.9])
        rho = 5.3
        
        x = stoichiometry.frac_mole_to_weight(n,MM)
        n2 = stoichiometry.frac_weight_to_mole(x,MM)
        np.testing.assert_allclose(n,n2)
        
        x = stoichiometry.frac_mole_to_volume(n,rho,MM)
        n2 = stoichiometry.frac_volume_to_mole(x,rho,MM)
        np.testing.assert_allclose(n,n2)
        
        x = stoichiometry.frac_weight_to_volume(n,rho)
        n2 = stoichiometry.frac_volume_to_weight(x,rho)
        np.testing.assert_allclose(n,n2)

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_stoichiometry("test_fraction"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

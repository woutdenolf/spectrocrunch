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

from .. import detectors
from .. import scintillators

import numpy as np


class test_objects(unittest.TestCase):

    def test_detectors(self):
        self.assertRaises(RuntimeError, detectors.AreaDetector.factory, "")
        o = detectors.AreaDetector.factory("pcoedge55")
        o = detectors.AreaDetector.factory("areadetector")

    def test_scintillators(self):
        self.assertRaises(RuntimeError, scintillators.Scintillator.factory, "", 0)
        o = scintillators.Scintillator.factory("GGG",10,{"Eu":0.03})

        #import matplotlib.pyplot as plt
        #plt.figure()
        #energy = np.linspace(2,9,100,dtype=float)
        #plt.plot(energy,o.transmission(energy))
        
        o = scintillators.Scintillator.factory("LSO",10,{"Tb":0.03})

        #plt.plot(energy,o.transmission(energy))
        #plt.show()

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_objects("test_detectors"))
    testSuite.addTest(test_objects("test_scintillators"))

    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

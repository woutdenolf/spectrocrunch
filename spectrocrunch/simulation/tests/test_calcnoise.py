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

from .. import calcnoise
from .. import materials

from ...materials.compoundfromformula import compoundfromformula as compound

import numpy as np
from uncertainties import unumpy

class test_calcnoise(unittest.TestCase):

    def test_ffnoise(self):
        I0 = 1e5
        energy = np.linspace(3,5,100)
        tframe = 0.07
        nframe = 100
        ndark = 30

        sample = materials.factory("Multilayer",material=compound("CaCO3",2.71),thickness=5)

        N,N0,D,D0 = calcnoise.id21_ffnoise(I0,energy,sample,\
                    tframe_data=tframe,nframe_data=nframe,\
                    tframe_flat=tframe,nframe_flat=nframe,\
                    nframe_dark=ndark)

        XAS = -unumpy.log((N-D)/(N0-D0))

        signal = unumpy.nominal_values(XAS)
        noise = unumpy.std_devs(XAS)

        #import matplotlib.pyplot as plt
        #plt.figure()
        #plt.plot(energy,noise/signal*100)
        #plt.show()
        
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_calcnoise("test_ffnoise"))

    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

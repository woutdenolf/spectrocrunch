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
from .. import lenses
from .. import materials
from ...materials.compoundfromformula import compound as compound

import numpy as np
from uncertainties import unumpy
from uncertainties import ufloat

class test_objects(unittest.TestCase):

    def _checkprop(self,o,**kwargs):
        # Both arrays
        energy = np.array([7.,7.])
        N = np.array([1e5,1e5])
        N = unumpy.uarray(N,np.sqrt(N))

        Nout = o.propagate(N,energy,**kwargs)

        self.assertEqual(len(Nout),2)
        tmp = unumpy.nominal_values(Nout)
        self.assertEqual(tmp[0],tmp[1])
        tmp = unumpy.std_devs(Nout)
        self.assertEqual(tmp[0],tmp[1])

        # N array
        energy = 7.
        N = np.array([1e5,1e5])
        N = unumpy.uarray(N,np.sqrt(N))

        Nout = o.propagate(N,energy,**kwargs)

        self.assertEqual(len(Nout),2)
        tmp = unumpy.nominal_values(Nout)
        self.assertEqual(tmp[0],tmp[1])
        tmp = unumpy.std_devs(Nout)
        self.assertEqual(tmp[0],tmp[1])

        # Energy array
        energy = np.array([7.,7.])
        N = 1e5
        N = ufloat(N,np.sqrt(N))

        Nout = o.propagate(N,energy,**kwargs)

        self.assertEqual(len(Nout),2)
        tmp = unumpy.nominal_values(Nout)
        self.assertEqual(tmp[0],tmp[1])
        tmp = unumpy.std_devs(Nout)
        self.assertEqual(tmp[0],tmp[1])

        # Not arrays
        energy = 7.
        N = 1e5
        N = ufloat(N,np.sqrt(N))
        
        Nout = o.propagate(N,energy,**kwargs)
        self.assertFalse(hasattr(Nout,"__iter__"))

    def test_detectors(self):
        self.assertRaises(RuntimeError, detectors.factory, "")

        o = detectors.factory("pcoedge55")
        o = detectors.pcoedge55()
        self._checkprop(o,tframe=2,nframe=10)

        o = detectors.factory("areadetector")
        self._checkprop(o,tframe=2,nframe=10)
    
    def test_lenses(self):
        self.assertRaises(RuntimeError, lenses.Lens.factory, "")

        o = lenses.factory("mitutoyoid21_10x")
        self._checkprop(o,nrefrac=1.1)

        o = lenses.factory("mitutoyoid21_20x")
        self._checkprop(o,nrefrac=1.1)

    def test_scintillators(self):
        self.assertRaises(RuntimeError, scintillators.factory, "", 0)

        o = scintillators.Scintillator.factory("GGG ID21",13)
        self._checkprop(o)

        #SNR = unumpy.nominal_values(N)/unumpy.std_devs(N)

        #import matplotlib.pyplot as plt
        #plt.figure()
        #energy = np.linspace(2,9,100,dtype=float)
        #plt.plot(energy,o.transmission(energy))
        
        o = scintillators.factory("LSO ID21",10)
        self._checkprop(o)

        #plt.plot(energy,o.transmission(energy))
        #plt.show()

    def test_materials(self):
        #self.assertRaises(RuntimeError, materials.Material.factory, "", 0)



        o = materials.factory("Ultralene",4)

        #for s in materials.Material.registry:
        #    print s
        
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_objects("test_detectors"))
    testSuite.addTest(test_objects("test_scintillators"))
    testSuite.addTest(test_objects("test_lenses"))
    testSuite.addTest(test_objects("test_materials"))

    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

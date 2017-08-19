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

from .. import areadetectors
from .. import scintillators
from .. import lenses
from .. import materials
from .. import diodes
from ...materials.compoundfromformula import compoundfromformula as compound
from ...materials.compoundfromname import compoundfromname as compoundname

import numpy as np
from uncertainties import unumpy
from uncertainties import ufloat

from ...resources import resource_filename
from scipy import constants
from scipy import interpolate
from ...math.linop import linop

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
        self.assertRaises(RuntimeError, areadetectors.factory, "noclassname")

        o = areadetectors.factory("PCO Edge 5.5")
        self._checkprop(o,tframe=2,nframe=10)
    
    def test_lenses(self):
        self.assertRaises(RuntimeError, lenses.factory, "noclassname")

        o = lenses.factory("Mitutoyo ID21 10x")
        self._checkprop(o,nrefrac=1.1)

        o = lenses.factory("Mitutoyo ID21 10x")
        self._checkprop(o,nrefrac=1.1)

    def test_scintillators(self):
        self.assertRaises(RuntimeError, scintillators.factory, "noclassname")
        for cname in scintillators.classes:
            self.assertRaises(RuntimeError, scintillators.factory, cname)

        o = scintillators.factory("GGG ID21",thickness=13)
        self._checkprop(o)

        #SNR = unumpy.nominal_values(N)/unumpy.std_devs(N)

        #import matplotlib.pyplot as plt
        #plt.figure()
        #energy = np.linspace(2,9,100,dtype=float)
        #plt.plot(energy,o.transmission(energy))
        
        o = scintillators.factory("LSO ID21",thickness=10)
        self._checkprop(o)

        #plt.plot(energy,o.transmission(energy))
        #plt.show()

    def test_materials(self):
        self.assertRaises(RuntimeError, materials.factory, "noclassname")
        o = materials.factory("multilayer",material=compoundname("ultralene"),thickness=4)
        self._checkprop(o)

    def _spec_calc_photons(self,ph_E,ph_I,ph_gain):
        """Photon calculation in spec

        Args:
            ph_E(array): keV
            ph_I(array): idet counts
            ph_gain(array): V/A
        """

        PH_DIST = 0
        ph_coeff = 0

        ph_I = ph_I * 10.**(-5-ph_gain)

        ptb = np.loadtxt(resource_filename('id21/ptb.dat'))
        fptb = interpolate.interp1d(ptb[:,0],ptb[:,1])
        ird = np.loadtxt(resource_filename('id21/ird.dat'))
        fird = interpolate.interp1d(ird[:,0],ird[:,1])

        ph_PTB = fptb(ph_E)
        ph_factor = fird(ph_E)

        ph_calib  = ph_factor * ph_PTB
        return ph_I / (ph_E * constants.elementary_charge * np.exp(-ph_coeff * PH_DIST) * ph_calib)

    def test_diodes(self):
        gain = 8
        I = np.arange(5,8)*1e5

        for energy in np.arange(3,9):
            for model in [True,False]:
                o = diodes.factory("sxmidet",model=model)
                o.setgain(10**gain)

                o2 = o.pndiode.op_cpstocurrent()*o.pndiode.op_currenttocps()
                self.assertEqual(o2.m,1.)
                self.assertEqual(o2.b,0.)

                o2 = o.pndiode.op_fluxtocurrent(energy)*o.pndiode.op_currenttoflux(energy)
                self.assertAlmostEqual(o2.m,1.)
                self.assertEqual(o2.b,0.)

                o2 = o.pndiode.op_fluxtocps(energy)*o.pndiode.op_cpstoflux(energy)
                self.assertAlmostEqual(o2.m,1.)
                self.assertEqual(o2.b,0.)

                np.testing.assert_array_almost_equal(o.fluxtocps(energy,o.cpstoflux(energy,I)),I)
                
                flux1 = self._spec_calc_photons(energy,I,gain)
                flux2 = o.cpstoflux(energy,I)
                if model:
                    for f1,f2 in zip(flux1,flux2):
                        np.testing.assert_approx_equal(f1,f2,significant=1)
                else:
                    np.testing.assert_array_almost_equal(flux1,flux2)

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_objects("test_detectors"))
    testSuite.addTest(test_objects("test_scintillators"))
    testSuite.addTest(test_objects("test_lenses"))
    testSuite.addTest(test_objects("test_materials"))
    testSuite.addTest(test_objects("test_diodes"))

    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

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

from .. import diode

from ...resources import resource_filename

from scipy import interpolate

from scipy import constants

from ... import ureg

import numpy as np

class test_objects(unittest.TestCase):                

    def _sxm_calc_photons(self,ph_E,ph_I,ph_gain):
        """Photon calculation in sxm

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

    def test_calibrateddiode(self):
        gain = 8
        I = np.arange(5,8)*ureg.Quantity(1e5,"1/s")
        
        o1 = diode.sxmidet(model=True)
        o2 = diode.sxmidet(model=False)
        self.assertAlmostEqual(o1._chargepersamplephoton(5.2).magnitude,o2._chargepersamplephoton(5.2).magnitude)

        for model in [True,False]:
            o = diode.sxmidet(model=model)
            o.setgain(ureg.Quantity(10**gain,'V/A'))

            o2 = o.op_cpstocurrent()*o.op_currenttocps()
            self.assertEqual(o2.m,1.)
            self.assertEqual(o2.b.magnitude,0.)

            for energy in np.arange(3,9):
                energy = ureg.Quantity(energy,"keV")
                o2 = o.op_fluxtocurrent(energy)*o.op_currenttoflux(energy)
                self.assertAlmostEqual(o2.m,1.)
                self.assertEqual(o2.b.magnitude,0.)

                o2 = o.op_fluxtocps(energy)*o.op_cpstoflux(energy)
                self.assertAlmostEqual(o2.m,1.)
                self.assertEqual(o2.b.magnitude,0.)

                np.testing.assert_array_almost_equal(o.fluxtocps(energy,o.cpstoflux(energy,I)),I)
                
                flux1 = self._sxm_calc_photons(energy.magnitude,I.magnitude,gain)
                flux2 = o.cpstoflux(energy,I)

                if model:
                    for f1,f2 in zip(flux1,flux2):
                        np.testing.assert_approx_equal(f1,f2.magnitude,significant=1)
                else:
                    np.testing.assert_allclose(flux1,flux2)

    def test_noncalibrateddiode(self):
        pass
    
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_objects("test_calibrateddiode"))
    testSuite.addTest(test_objects("test_noncalibrateddiode"))
    
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

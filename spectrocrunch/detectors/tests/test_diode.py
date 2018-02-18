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

class test_diode(unittest.TestCase):                

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
        
        o1 = diode.SXM_IDET(model=True)
        o2 = diode.SXM_IDET(model=False)
        self.assertAlmostEqual(o1._chargepersamplephoton(5.2).magnitude,o2._chargepersamplephoton(5.2).magnitude)

        for model in [True,False]:
            o = diode.SXM_IDET(model=model)
            o.setgain(ureg.Quantity(10**gain,'V/A'))

            o2 = o.op_cpstocurrent()*o.op_currenttocps()
            self.assertEqual(o2.m,1.)
            self.assertEqual(o2.b.magnitude,0.)
            self.assertEqual(o2.m.units,ureg.dimensionless)
            self.assertEqual(o2.b.units,ureg.ampere)

            for energy in np.arange(3,9):
                o2 = o.op_fluxtocurrent(energy)*o.op_currenttoflux(energy)
                self.assertAlmostEqual(o2.m,1.)
                self.assertAlmostEqual(o2.b.magnitude,0.)
                self.assertEqual(o2.m.units,ureg.dimensionless)
                self.assertEqual(o2.b.units,ureg.ampere)

                o2 = o.op_fluxtocps(energy)*o.op_cpstoflux(energy)
                self.assertAlmostEqual(o2.m,1.)
                self.assertAlmostEqual(o2.b.magnitude,0.)
                self.assertEqual(o2.m.units,ureg.dimensionless)
                self.assertEqual(o2.b.units,ureg.hertz)
                
                np.testing.assert_array_almost_equal(o.fluxtocps(energy,o.cpstoflux(energy,I)),I)
                
                flux1 = self._sxm_calc_photons(energy,I.magnitude,gain)
                flux2 = o.cpstoflux(energy,I)

                if model:
                    for f1,f2 in zip(flux1,flux2):
                        np.testing.assert_approx_equal(f1,f2.magnitude,significant=1)
                else:
                    np.testing.assert_allclose(flux1,flux2,rtol=0.03) # 3% difference with spec

    def test_noncalibrateddiode(self):
    
        for caliboption in ["solidangle","thickness","optics"]:
            o = diode.SXM_IODET1(optics=True)
            o.thickness = 1e-4
            o.geometry.solidangle = 0.1
            o.optics.reset_transmission()

            cps = np.linspace(1e4,1e5,10)
            
            energy = 7.
            energy2 = 7.2
            energy3 = 7.1
            sampleflux = o.cpstoflux(energy,cps)
            sampleflux2 = o.cpstoflux(energy2,cps)
            sampleflux3 = o.cpstoflux(energy3,cps)

            if caliboption=="solidangle":
                # Same counts but the flux is m times higher -> diode yield (solid angle) is m times lower
                m = 2.
                sa = o.geometry.solidangle
                o.calibrate(cps,sampleflux*m,energy,caliboption=caliboption)
                np.testing.assert_allclose(sampleflux*m,o.cpstoflux(energy,cps))
                np.testing.assert_allclose(o.geometry.solidangle,sa/m)
            elif caliboption=="optics":
                # Same counts but the flux is m times lower -> optics transmission is m times lower
                m1 = 0.5
                m2 = 0.25
                m3 = (m1+m2)/2.

                o.calibrate(cps,sampleflux*m1,energy,caliboption=caliboption)
                o.calibrate(cps,sampleflux2*m2,energy2,caliboption=caliboption)

                np.testing.assert_allclose(o.optics.transmission(energy),m1)
                np.testing.assert_allclose(o.optics.transmission(energy2),m2)
                np.testing.assert_allclose(o.optics.transmission(energy3),m3)
                
                np.testing.assert_allclose(sampleflux*m1,o.cpstoflux(energy,cps))
                np.testing.assert_allclose(sampleflux2*m2,o.cpstoflux(energy2,cps))
                np.testing.assert_allclose(sampleflux3*m3,o.cpstoflux(energy3,cps))
            elif caliboption=="thickness":
                m = 2.
                o.calibrate(cps,sampleflux*m,energy,caliboption=caliboption)
                np.testing.assert_allclose(sampleflux*m,o.cpstoflux(energy,cps))
            
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_diode("test_calibrateddiode"))
    testSuite.addTest(test_diode("test_noncalibrateddiode"))
    
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
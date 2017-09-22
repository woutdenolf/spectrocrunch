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

from ..element import element
from ..element import shell
from ..element import fluoline

import numpy as np
import xraylib

class test_element(unittest.TestCase):

    def test_fluoline(self):
        self.assertEqual(fluoline.factory(None,"K"),None)
        self.assertEqual(fluoline.factory([],"K"),[])
        self.assertTrue(fluoline.factory(xraylib.KA_LINE,"L") is None)
        self.assertEqual(fluoline.factory(xraylib.KA1_LINE,"K"),fluoline.factory([xraylib.KA1_LINE],"K"))
        self.assertEqual(fluoline.factory(xraylib.KA_LINE,"K"),fluoline.factory([xraylib.KA1_LINE,xraylib.KA2_LINE,xraylib.KA3_LINE],"K"))
        self.assertEqual(fluoline.factory(xraylib.KA_LINE,"K"),fluoline.factory([xraylib.KA_LINE,xraylib.LA_LINE],"K"))
        self.assertNotEqual(fluoline.factory(xraylib.KA_LINE,"K"),fluoline.factory([xraylib.KA1_LINE,xraylib.KA3_LINE],"K"))

    def test_shell(self):
        self.assertEqual(shell(xraylib.K_SHELL).radrate(26),[])
        self.assertEqual(shell(xraylib.K_SHELL,fluolines=[]).radrate(26),[1])
        self.assertEqual(len(shell(xraylib.K_SHELL,fluolines=xraylib.KA_LINE).radrate(26)),3)
        self.assertEqual(len(shell(xraylib.K_SHELL,fluolines=xraylib.LA_LINE).radrate(26)),0)
        self.assertEqual(len(shell(xraylib.K_SHELL,fluolines=[xraylib.LA_LINE,xraylib.KA_LINE]).radrate(26)),3)
        self.assertEqual(sum(shell(xraylib.K_SHELL,fluolines=[xraylib.KA_LINE,xraylib.KB_LINE]).radrate(26)),1)

    def test_fluo(self):
        e = element("Fe")

        mu1 = 0.
        mu2 = 0.
        for s in shell._all:
            mu1 += xraylib.CS_Photo_Partial(26,s,8.)
            mu2 += xraylib.CS_Photo_Partial(26,s,8.)*xraylib.FluorYield(26,s)
        e.markasabsorber("Fe",shells=shell._all.keys(),fluolines=[])
        mu3 = e.partial_mass_abs_coeff(8.)
        mu4 = e.fluorescence_cross_section(8.)
        self.assertEqual(mu1,mu3)
        self.assertEqual(mu2,mu4)

        mu1 = xraylib.CS_Photo_Partial(26,xraylib.K_SHELL,8.)
        mu2 = xraylib.CS_Photo_Partial(26,xraylib.K_SHELL,8.)*xraylib.FluorYield(26,xraylib.K_SHELL)
        e.markasabsorber("Fe",shells=xraylib.K_SHELL,fluolines=[])
        mu3 = e.partial_mass_abs_coeff(8.)
        mu4 = e.fluorescence_cross_section(8.)
        self.assertEqual(mu1,mu3)
        self.assertEqual(mu2,mu4)
        e.markasabsorber("Fe",shells=xraylib.K_SHELL,fluolines=[xraylib.KA_LINE,xraylib.KB_LINE,xraylib.KA_LINE,xraylib.KB_LINE])
        mu3 = e.partial_mass_abs_coeff(8.)
        mu4 = e.fluorescence_cross_section(8.)
        self.assertEqual(mu1,mu3)
        self.assertEqual(mu2,mu4)

        mu1 = xraylib.CS_Photo_Partial(26,xraylib.K_SHELL,8.)
        mu2 = xraylib.CS_Photo_Partial(26,xraylib.K_SHELL,8.)*xraylib.FluorYield(26,xraylib.K_SHELL)*xraylib.RadRate(26,xraylib.KA_LINE)
        e.markasabsorber("Fe",shells=xraylib.K_SHELL,fluolines=[xraylib.KA1_LINE,xraylib.KA2_LINE,xraylib.KA3_LINE,xraylib.LA_LINE])
        mu3 = e.partial_mass_abs_coeff(8.)
        mu4 = e.fluorescence_cross_section(8.)
        self.assertEqual(mu1,mu3)
        self.assertEqual(mu2,mu4)

        mu1 = e.mass_att_coeff(8.)
        mu2 = e.mass_abs_coeff(8.) + e.scattering_cross_section(8.)
        self.assertEqual(mu1,mu2)
        mu2 = e.mass_abs_coeff(8.) + e.rayleigh_cross_section(8.) + e.compton_cross_section(8.)
        self.assertEqual(mu1,mu2)

        e.markasabsorber("Fe",shells=shell._all.keys(),fluolines=[])
        energy = np.linspace(5,10,500)
        mu1 = e.partial_mass_abs_coeff(energy)
        mu2 = e.mass_abs_coeff(energy)
        #np.testing.assert_almost_equal(mu1,mu2) # Auger is negligible

        #import matplotlib.pyplot as plt
        #plt.plot(energy,mu1)
        #plt.plot(energy,mu2)
        #plt.show()

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_element("test_fluoline"))
    testSuite.addTest(test_element("test_shell"))
    testSuite.addTest(test_element("test_fluo"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

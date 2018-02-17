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

from .. import element
from .. import xrayspectrum
from ... import xraylib
from ... import ureg

import numpy as np
from scipy import integrate

class test_element(unittest.TestCase):

    def _linesequal(self,lines1,lines2):
        self.assertEqual(sorted(lines1),sorted(lines2))

    def test_fluoline(self):
    
        #print [xrayspectrum.FluoLine.getlinename(code) for code in  sorted([xrayspectrum.FluoLine.getlinecode(name) for name in xrayspectrum.FluoLine.decompose("L2P23")])]
        
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="K")),29-2)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="L1")),58-29-3)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="L2")),85-58-1)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="L3")),113-85-3)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="M1")),136-113)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="M2")),158-136)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="M3")),180-158)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="M4")),200-180)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="M5")),219-200)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="N1")),237-219)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="N2")),254-237)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="N3")),270-254)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="N4")),285-270)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="N5")),299-285)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="N6")),312-299)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="N7")),324-312)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="O1")),335-324)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="O2")),345-335)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="O3")),354-345)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="O4")),362-354)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="O5")),369-362)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="O6")),372-369)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="O7")),374-372)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="P1")),378-374)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="P2")),381-378)
        self.assertEqual(len(xrayspectrum.FluoLine.factory(shells="P3")),383-381)
        
        energybounds = None
        fluolines = []
        lines = list(set(range(-29,0)) - set([xraylib.KO_LINE,xraylib.KP_LINE]))
        linestest = [xrayspectrum.FluoLine(line) for line in lines]
        self._linesequal(xrayspectrum.FluoLine.factory(shells="K",fluolines=fluolines,energybounds=energybounds),linestest)
        
        energybounds = None
        lines = [xraylib.KA1_LINE,xraylib.KA2_LINE,xraylib.KA3_LINE,xraylib.KA_LINE]
        fluolines = lines + [xraylib.LA_LINE]
        linestest = [xrayspectrum.FluoLine(line) for line in lines]
        self._linesequal(xrayspectrum.FluoLine.factory(shells="K",fluolines=fluolines,energybounds=energybounds),linestest)

        linestest = []
        self._linesequal(xrayspectrum.FluoLine.factory(shells="L1",fluolines=fluolines,energybounds=energybounds),linestest)
        
        linestest = [xrayspectrum.FluoLine(xraylib.LA_LINE)]
        self._linesequal(xrayspectrum.FluoLine.factory(shells="L3",fluolines=fluolines,energybounds=energybounds),linestest)
        
        energybounds = (element.elementZ('Fe'),0,6.395)
        lines = [xraylib.KA2_LINE,xraylib.KA3_LINE]
        linestest = [xrayspectrum.FluoLine(line) for line in lines]
        self._linesequal(xrayspectrum.FluoLine.factory(shells="K",fluolines=fluolines,energybounds=energybounds),linestest)
        
    def test_shell(self):
        Z = 50
        for shell in ['K','L1','L2','L3','M1','M2','M3','M4','M5']:
            shell = xrayspectrum.Shell(shell,fluolines=[])
            np.testing.assert_allclose(sum(shell.radrate(Z)),1,rtol=1e-5)
    
        self.assertEqual(xrayspectrum.Shell(xraylib.K_SHELL,fluolines=None).radrate(26),[1])
        self.assertEqual(len(xrayspectrum.Shell(xraylib.K_SHELL,fluolines=[]).radrate(26)),27)
        self.assertEqual(sum(xrayspectrum.Shell(xraylib.K_SHELL,fluolines=[]).radrate(26)),1)
        self.assertEqual(sum(xrayspectrum.Shell(xraylib.K_SHELL,fluolines=[xraylib.KA_LINE,xraylib.KB_LINE]).radrate(26)),1)
        self.assertEqual(len(xrayspectrum.Shell(xraylib.K_SHELL,fluolines=[xraylib.KA1_LINE,xraylib.LA_LINE,xraylib.KA_LINE]).radrate(26)),2)        
        self.assertEqual(len(xrayspectrum.Shell(xraylib.K_SHELL,fluolines=[xraylib.LA_LINE]).radrate(26)),0)
        self.assertEqual(len(xrayspectrum.Shell(xraylib.L3_SHELL,fluolines=[xraylib.LA1_LINE,xraylib.KA_LINE]).radrate(26)),1)

    def test_fluo(self):
        e = element.Element("Sn")

        for s in xraylib.code_to_shell:
            mu1 = xraylib.CS_Photo_Partial(50,s,8.)
            mu2 = xraylib.CS_Photo_Partial(50,s,8.)*xraylib.FluorYield(50,s)
            
            e.markabsorber("Sn",shells=[s],fluolines=None) # all lines implicitely
            mu3 = e.partial_mass_abs_coeff(8.)
            mu4 = e.fluorescence_cross_section(8.)
            self.assertEqual(mu1,mu3)
            self.assertEqual(mu2,mu4)
            
            e.markabsorber("Sn",shells=[s],fluolines=[]) # all lines explicitely
            mu3 = e.partial_mass_abs_coeff(8.)
            mu4 = e.fluorescence_cross_section(8.)
            self.assertEqual(mu1,mu3)
            np.testing.assert_allclose(mu2,mu4,rtol=1e-5,atol=2e-5)

        mu1 = 0.
        mu2 = 0.
        for s in xraylib.code_to_shell:
            mu1 += xraylib.CS_Photo_Partial(50,s,8.)
            mu2 += xraylib.CS_Photo_Partial(50,s,8.)*xraylib.FluorYield(50,s)
            
        e.markabsorber("Sn",shells=None,fluolines=None) # all lines implicitely
        mu3 = e.partial_mass_abs_coeff(8.)
        mu4 = e.fluorescence_cross_section(8.)
        self.assertEqual(mu1,mu3)
        np.testing.assert_allclose(mu2,mu4,rtol=1e-5,atol=2e-5)
        
        e.markabsorber("Sn",shells=None,fluolines=[]) # all lines explicitely
        mu3 = e.partial_mass_abs_coeff(8.)
        mu4 = e.fluorescence_cross_section(8.)
        self.assertEqual(mu1,mu3)
        np.testing.assert_allclose(mu2,mu4,rtol=1e-5,atol=2e-5)
            
        e = element.Element("Fe")
        
        mu1 = xraylib.CS_Photo_Partial(26,xraylib.K_SHELL,8.)
        mu2 = xraylib.CS_Photo_Partial(26,xraylib.K_SHELL,8.)*xraylib.FluorYield(26,xraylib.K_SHELL)
        e.markabsorber("Fe",shells=xraylib.K_SHELL,fluolines=[])
        mu3 = e.partial_mass_abs_coeff(8.)
        mu4 = e.fluorescence_cross_section(8.)
        self.assertEqual(mu1,mu3)
        self.assertEqual(mu2,mu4)
        e.markabsorber("Fe",shells=xraylib.K_SHELL,fluolines=[xraylib.KA_LINE,xraylib.KB_LINE])
        mu3 = e.partial_mass_abs_coeff(8.)
        mu4 = e.fluorescence_cross_section(8.)
        self.assertEqual(mu1,mu3)
        self.assertEqual(mu2,mu4)

        mu1 = xraylib.CS_Photo_Partial(26,xraylib.K_SHELL,8.)
        mu2 = xraylib.CS_Photo_Partial(26,xraylib.K_SHELL,8.)*xraylib.FluorYield(26,xraylib.K_SHELL)*xraylib.RadRate(26,xraylib.KA_LINE)
        e.markabsorber("Fe",shells=xraylib.K_SHELL,fluolines=[xraylib.KA1_LINE,xraylib.KA2_LINE,xraylib.KA3_LINE,xraylib.LA_LINE])
        mu3 = e.partial_mass_abs_coeff(8.)
        mu4 = e.fluorescence_cross_section(8.)
        self.assertEqual(mu1,mu3)
        self.assertEqual(mu2,mu4)

        mu1 = e.mass_att_coeff(8.)
        mu2 = e.mass_abs_coeff(8.) + e.scattering_cross_section(8.)
        self.assertEqual(mu1,mu2)
        mu2 = e.mass_abs_coeff(8.) + e.rayleigh_cross_section(8.) + e.compton_cross_section(8.)
        self.assertEqual(mu1,mu2)

        e.markabsorber("Fe",shells=None,fluolines=[])
        energy = np.linspace(5,10,500)
        mu1 = e.partial_mass_abs_coeff(energy)
        mu2 = e.mass_abs_coeff(energy)
        np.testing.assert_almost_equal(mu1,mu2) # Auger is negligible

        #import matplotlib.pyplot as plt
        #plt.plot(energy,mu1)
        #plt.plot(energy,mu2)
        #plt.show()
        
    def test_comparable(self):
        self.assertEqual(element.Element("Ca"),element.Element("Ca"))
        self.assertNotEqual(element.Element("Ca"),element.Element("C"))
        self.assertEqual(element.Element("Ca"),"Ca")
        self.assertNotEqual(element.Element("Ca"),"C")
    
    def test_diffcs_elastic(self):
        e = element.Element("Fe")
        
        energy = 30
        theta = np.radians(45)
        np.testing.assert_allclose(xraylib.FF_Rayl(e.Z,xraylib.MomentTransf(energy,theta)),e.scatfact_classic_re(energy,theta=theta),rtol=1e-5)
        
        # Elastic scattering cross-section (muR) as an integral of the atomic form factor (f)
        #
        #from sage.all import *
        #from sage.symbolic.integration.integral import definite_integral
        #
        #var('P')
        #var('theta')
        #var('phi')
        #
        #K(phi) = (1-P*cos(2*phi))/2
        #ElasticDiff2(theta,phi) = 1-K(phi)+K(phi)*(cos(theta))**2
        #
        #ElasticDiff1(theta) = definite_integral(ElasticDiff2(theta,phi),phi,0,2*pi)
        #
        #print ElasticDiff1(theta)
        #
        #assert(definite_integral(ElasticDiff1(theta)*sin(theta),theta,0,pi)==8*pi/3)
        #
        # muR(energy) = re^2 . NA/MM . pi . int_0^pi [(1+cos(theta)^2) . f(energy,theta)^2 . sin(theta) . dtheta]
        
        energy = np.linspace(4,30,2)
        fim2 = e.scatfact_im(energy)**2
        
        def integrand1(theta,energy,a):
            return (1+np.cos(theta)**2)*np.sin(theta)*(e.scatfact_re(energy,theta=theta)**2+a)
            
        def integrand2(theta,energy):
            return (1+np.cos(theta)**2)*np.sin(theta)*e.scatfact_classic_re(energy,theta=theta)**2
        
        fac = (np.pi*ureg.re**2*ureg.avogadro_number/ureg.Quantity(e.MM,'g/mol')).to("cm^2/g").magnitude
        
        #cs1 = [fac*integrate.quad(integrand1, 0, np.pi, args = (energy[i],fim2[i]))[0]  for i in range(len(energy))]
        
        cs2 = [fac*integrate.quad(integrand2, 0, np.pi, args = (E))[0]  for E in energy]
        
        cs3 = e.rayleigh_cross_section(energy)

        np.testing.assert_allclose(cs2,cs3,rtol=1e-1)
        
        #import matplotlib.pyplot as plt
        #plt.plot(energy,cs1,label='Int Full')
        #plt.plot(energy,cs2,label='Int Classic')
        #plt.plot(energy,cs3,label='Tab')
        #plt.legend()
        #plt.show()
            
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_element("test_fluoline"))
    testSuite.addTest(test_element("test_shell"))
    testSuite.addTest(test_element("test_fluo"))
    testSuite.addTest(test_element("test_comparable"))
    testSuite.addTest(test_element("test_diffcs_elastic"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

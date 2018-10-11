# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
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
import cmath

from .. import polarization
from ...utils import instance

import numpy as np
from scipy import integrate

class test_polarization(unittest.TestCase):

    def _equal_params(self,params1,params2):
        for k,v in params1.items():
            if instance.isstring(v):
                self.assertEqual(v,params2[k])
            else:
                np.testing.assert_allclose(v,params2[k])
                
    def _gen_jones(self,n=20):
        x = np.random.uniform(low=-10,high=10,size=4*n).reshape((n,4))
        for xi in x:
            yield polarization.Jones(xi[0]+xi[1]*1j,xi[2]+xi[3]*1j)

    def _gen_stokes(self,n=20):
        x = np.random.uniform(low=-10,high=10,size=3*n).reshape((n,3))
        for xi in x:
            S0 = np.sqrt(sum(xi[1:]**2))*np.random.uniform(low=1,high=1.5)
            yield polarization.Stokes(S0,*xi)
    
    def test_convert_representation(self):
        f1 = lambda x,attr:getattr(x,attr)
        f2 = lambda x,attr:getattr(x,attr)%360
        attrs = {"coherency_matrix":f1,"dop":f1,"dolp":f1,"docp":f1,"hdolp":f1,"polangle":f2}
        
        for J1 in self._gen_jones():
            S1 = J1.to_stokes()
            J2 = S1.to_jones()
            S2 = J2.to_stokes()
            J3 = S2.to_jones()

            self._equal_params(J2.to_params(),J3.to_params())
            self._equal_params(S1.to_params(),S2.to_params())

            self.assertEqual(J1.dop,1)
            
            for attr,f in attrs.items():
                a = f(J1,attr)
                np.testing.assert_allclose(a,f(S1,attr))
                np.testing.assert_allclose(a,f(J2,attr))
                np.testing.assert_allclose(a,f(S2,attr))
                np.testing.assert_allclose(a,f(J3,attr))

            np.testing.assert_allclose(J1.norm,J2.norm)
            np.testing.assert_allclose(J1.phase_difference%360,J2.phase_difference%360)
            np.testing.assert_allclose(J2.to_numpy(),J3.to_numpy())
            np.testing.assert_allclose(S1.to_numpy(),S2.to_numpy())
            np.testing.assert_allclose(S1.to_numpy(),S2.to_numpy())

    def test_stokes(self):
        for S in self._gen_stokes():
            tmp = S.decompose()
            Spol,Sunpol = tmp["pol"],tmp["unpol"]
            np.testing.assert_allclose(S.intensity,S.intensity_polarized+S.intensity_unpolarized)
            np.testing.assert_allclose(S.intensity_polarized,Spol.intensity)
            np.testing.assert_allclose(S.intensity_unpolarized,Sunpol.intensity)
            np.testing.assert_allclose(S.dop,S.intensity_polarized/S.intensity)
            
            np.testing.assert_allclose(S.coherency_matrix,Spol.coherency_matrix+Sunpol.coherency_matrix)
            
            J = S.to_jones(allowloss=True)
            np.testing.assert_allclose(J.intensity,Spol.intensity)
            
            S2 = polarization.Stokes.from_params(**S.to_params())
            np.testing.assert_allclose(S.to_numpy(),S2.to_numpy())
    
    def test_jones(self):
        for J in self._gen_jones():
            np.testing.assert_allclose(J.to_numpy(),J.to_stokes().to_jones(phase0=J.phase0).to_numpy())
            np.testing.assert_allclose(J.coherency_matrix.trace(),J.norm**2)
            
            J2 = polarization.Jones.from_params(**J.to_params())
            np.testing.assert_allclose(J.to_numpy(),J2.to_numpy())
            
            J.plot_efield(animate=True)

    def test_intensity(self):
        for J in self._gen_jones():
            S = J.to_stokes()
            Jparams = J.to_params()
            Sparams = S.to_params()
            IJ,IS = np.random.uniform(low=1,high=10,size=2)
            J.intensity = IJ
            S.intensity = IS
            Jparams["intensity"] = IJ
            Sparams["intensity"] = IS
            self._equal_params(J.to_params(),Jparams)
            self._equal_params(S.to_params(),Sparams)

        for S in self._gen_stokes():
            Sparams = S.to_params()
            IS = np.random.uniform(low=1,high=10)
            S.intensity = IS
            Sparams["intensity"] = IS
            self._equal_params(S.to_params(),Sparams)

    def test_rotate(self):
        for J1 in self._gen_jones():
            S1 = J1.to_stokes()
            
            azimuth = np.random.uniform(low=0,high=2*np.pi) # change-of-frame
            
            J2 = J1.rotate(azimuth)
            S2 = S1.rotate(azimuth)
            
            self._equal_params(S2.to_params(),J2.to_stokes().to_params())
            
            R = polarization.JonesMatrixRotation(-azimuth)
            Ri = polarization.JonesMatrixRotation(azimuth)
            
            np.testing.assert_allclose(R.dot(J1.coherency_matrix).dot(Ri),J2.coherency_matrix)
            np.testing.assert_allclose(R.dot(S1.coherency_matrix).dot(Ri),S2.coherency_matrix)
            
    def test_thomson(self):
        for J1 in self._gen_jones():
            S1 = J1.to_stokes()
            
            azimuth = np.random.uniform(low=0,high=2*np.pi)
            polar = np.random.uniform(low=0,high=np.pi)
            
            J2 = J1.thomson_scattering(azimuth,polar)
            S2 = S1.thomson_scattering(azimuth,polar)

            self._equal_params(S2.to_params(),J2.to_stokes().to_params())
            
            angle = polarization.ThomsonRotationAngle(azimuth) # change-of-frame
            R = polarization.JonesMatrixRotation(-angle)
            Ri = polarization.JonesMatrixRotation(angle)
            Mth = polarization.JonesMatrixThomson(polar)
            Mthi = Mth

            np.testing.assert_allclose(Mth.dot(R).dot(J1.coherency_matrix).dot(Ri).dot(Mthi),J2.coherency_matrix)
            np.testing.assert_allclose(Mth.dot(R).dot(S1.coherency_matrix).dot(Ri).dot(Mthi),S2.coherency_matrix)
            np.testing.assert_allclose(S2.intensity,S1.thomson_intensity(azimuth,polar))
            
            integrand = lambda azimuth,polar: S1.thomson_intensity(np.degrees(azimuth),np.degrees(polar))*np.sin(polar)
            thomsonsc = integrate.dblquad(integrand, 0, np.pi, lambda x:0, lambda x:2*np.pi)[0]/S1.intensity
            np.testing.assert_allclose(thomsonsc,8*np.pi/3)

    def test_compton(self):
        for S1 in self._gen_stokes():
            
            azimuth = np.random.uniform(low=0,high=2*np.pi)
            polar = np.random.uniform(low=0,high=np.pi)
            energy = np.random.uniform(low=5.,high=20.)
            
            S2 = S1.compton_scattering(azimuth,polar,energy)
            
            np.testing.assert_allclose(S2.intensity,S1.compton_intensity(azimuth,polar,energy))
            
            
def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_polarization("test_jones"))
    testSuite.addTest(test_polarization("test_stokes"))
    testSuite.addTest(test_polarization("test_convert_representation"))
    testSuite.addTest(test_polarization("test_intensity"))
    testSuite.addTest(test_polarization("test_rotate"))
    testSuite.addTest(test_polarization("test_thomson"))
    testSuite.addTest(test_polarization("test_compton"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

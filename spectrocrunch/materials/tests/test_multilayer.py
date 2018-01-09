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

from .. import multilayer
from .. import compoundfromformula
from .. import compoundfromname
from .. import mixture
from .. import types
from .. import utils
from .. import xrayspectrum
from ...geometries import xrf as xrfgeometries
from ...detectors import xrf as xrfdetectors

import numpy as np
import scipy.integrate
import scipy.special
import xraylib

class test_multilayer(unittest.TestCase):

    def test_expi(self):
        for x in [1e-4,1,10,np.inf]: # x > 0
            int1 = scipy.integrate.quad(lambda a: np.exp(-x/np.cos(a))*np.tan(a), 0, np.pi/2)
            int2 = -scipy.special.expi(-x)
            int3 = scipy.special.exp1(x)
            self.assertGreaterEqual(int2,int1[0]-int1[1])
            self.assertLessEqual(int2,int1[0]+int1[1])
            self.assertGreaterEqual(int3,int1[0]-int1[1])
            self.assertLessEqual(int3,int1[0]+int1[1])
            
    def _multilayer(self):
        geometry = xrfgeometries.factory("Geometry",detectorposition=0, anglein=62, angleout=49)
        detector = xrfdetectors.factory("Detector",activearea=0.50, geometry=geometry)
        detector.solidangle = 4*np.pi*0.1
    
        c1 = compoundfromformula.CompoundFromFormula("CaCO3",2.71,name="calcite")
        c2 = compoundfromformula.CompoundFromFormula("Fe2O3",5.3,name="hematite")
        c3 = compoundfromformula.CompoundFromFormula("PbCO3",6.53,name="cerussite")

        l = [c1,mixture.Mixture([c2,c3],[1,1],types.fractionType.mole)]
        thickness = [10,20]
        o = multilayer.Multilayer(material=l, thickness=thickness, detector=detector)
        
        return o,thickness

    def test_str(self):
        o,thickness = self._multilayer()
        s = str(o).split("\n")
        for i,layer in enumerate(o):
            self.assertTrue(str(layer) in s[i+1])
            self.assertTrue(str(thickness[i]) in s[i+1])
            
    def test_transmission(self):

        geometry = xrfgeometries.factory("geometry",anglein=90,angleout=45,detectorposition=0)
        detector = xrfdetectors.factory("leia",geometry=geometry)
    
        o = multilayer.Multilayer(\
                material=[compoundfromname.compoundfromname("hematite"),compoundfromname.compoundfromname("hydrocerussite"),compoundfromname.compoundfromname("calcite")],\
                thickness=[9e-4,9e-4,9e-4],detector=detector)
    
        with o.cachectx("layerinfo"):
        
            for energy in [8,np.array([8,8.5])]:
            
                with o.cachectx("attenuationinfo",energy):
                    nenergy = np.asarray(energy).size
                    
                    rho = o.density[:,np.newaxis]
                    d = o.xraythickness[:,np.newaxis]
                    mu = o.mass_att_coeff(energy)
                    murhod = mu*rho*d

                    nsub = 3
                    n = nsub*o.nlayers+1
                    A = np.zeros((n,nenergy))
                    np.cumsum(np.repeat(murhod,nsub,axis=0)/nsub,out=A[1:,:],axis=0)
                    
                    z = np.zeros(n)
                    np.cumsum(np.repeat(d,nsub)/nsub,out=z[1:])

                    np.testing.assert_allclose(A,o._cum_attenuation(z,energy))
                    
                    for i in range(n):
                        np.testing.assert_allclose(A[i,:],o._cum_attenuation(z[i],energy)[0,:])
                        cosaij = np.ones(n)
                        cosaij[0:i] = -1
                        
                        T1 = np.exp(-np.abs(A-A[i,np.newaxis,:]))
                        T2 = o._transmission(z[i],z,cosaij,energy)
                        np.testing.assert_allclose(T1,T2)
                        
                        for j in range(n):
                            T1 = np.exp(-abs(A[j,:]-A[i,:]))
                            T2 = o._transmission(z[i],z[j],cosaij[j],energy)[0,:]
                            np.testing.assert_allclose(T1,T2)
    
    def test_primary_fluorescence(self):
        c1 = compoundfromname.compoundfromname("hydrocerussite")
        c2 = compoundfromname.compoundfromname("hematite")
        c3 = compoundfromname.compoundfromname("calcite")
        
        geometry = xrfgeometries.factory("sdd120",detectorposition=-15.)
        detector = xrfdetectors.factory("leia",geometry=geometry)
        
        o = multilayer.Multilayer(\
                    material=[c1,c2,mixture.Mixture([c2,c3],[0.5,0.5],types.fractionType.weight),c3],\
                    thickness=[9e-4,12e-4,5e-4,15e-4],detector=detector)
        
        o.markabsorber("Fe",shells=xraylib.K_SHELL,fluolines=xrayspectrum.FluoLine.factory(energybounds = ("Fe",2,8)))
        o.markabsorber("Ca",shells=xraylib.K_SHELL,fluolines=xrayspectrum.FluoLine.factory(energybounds = ("Ca",2,8)))
        
        energy = np.array([8,8.5])
        #o._fluorescence(energy)
        
    def test_transmission2(self):
        o,thickness = self._multilayer()
        
        energy = np.linspace(7,8,5)

        T1 = o.transmission(energy)
        T2 = o.transmission(energy,decomposed=True)
        
        T3 = np.ones_like(energy)
        for mat in T2:
            for cs in utils.elemental_crosssections(mat["cs"]).values():
                T3 *= np.exp(-mat["density"]*mat["thickness"]*cs)

        np.testing.assert_allclose(T1,T3)
        
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_multilayer("test_expi"))
    testSuite.addTest(test_multilayer("test_str"))
    testSuite.addTest(test_multilayer("test_transmission"))
    testSuite.addTest(test_multilayer("test_transmission2"))
    testSuite.addTest(test_multilayer("test_primary_fluorescence"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

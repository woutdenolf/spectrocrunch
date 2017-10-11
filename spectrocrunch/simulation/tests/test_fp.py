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

import scipy.integrate as integrate
import scipy.special
import numpy as np
import xraylib

from .. import materials
from ...materials.compoundfromname import compoundfromname as compoundname
from ...materials.mixture import Mixture as mixture
from ...materials import element
from ...materials.types import fractionType

class test_fp(unittest.TestCase):

    def test_expi(self):
        for x in [1e-4,1,10,np.inf]: # x > 0
            int1 = integrate.quad(lambda a: np.exp(-x/np.cos(a))*np.tan(a), 0, np.pi/2)
            int2 = -scipy.special.expi(-x)
            int3 = scipy.special.exp1(x)
            self.assertGreaterEqual(int2,int1[0]-int1[1])
            self.assertLessEqual(int2,int1[0]+int1[1])
            self.assertGreaterEqual(int3,int1[0]-int1[1])
            self.assertLessEqual(int3,int1[0]+int1[1])
            
    def test_transmission(self):

        o = materials.factory("multilayer",\
                material=[compoundname("hematite"),compoundname("hydrocerussite"),compoundname("calcite")],\
                thickness=[9,12,15],\
                anglein=90-62,\
                angleout=90+49)
    
        for energy in [8,np.array([8,8.5])]:
        
            nenergy = np.asarray(energy).size
            
            rho = o.layerproperties("density").reshape((o.nlayers,1))
            d = o.thickness.reshape((o.nlayers,1))
            mu = o.layercs("mass_att_coeff",energy)
            mu = mu.reshape((o.nlayers,nenergy))  
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
                
                T1 = np.exp(-  abs(A-A[i,:].reshape((1,nenergy))) )
                T2 = o._transmission(z[i],z,cosaij,energy)
                np.testing.assert_allclose(T1,T2)
                
                for j in range(n):
                    T1 = np.exp(-abs(A[j,:]-A[i,:]))
                    T2 = o._transmission(z[i],z[j],cosaij[j],energy)[0,:]
                    np.testing.assert_allclose(T1,T2)
    
    def test_primary_fluorescence(self):
        c1 = compoundname("hydrocerussite")
        c2 = compoundname("hematite")
        c3 = compoundname("calcite")
        
        o = materials.factory("multilayer",\
                    material=[c1,c2,mixture([c2,c3],[0.5,0.5],fractionType.weight),c3],\
                    thickness=[9,12,5,15],\
                    anglein=90-62,\
                    angleout=90+49)
        
        
        o.markabsorber("Fe",shells=xraylib.K_SHELL,fluolines=element.FluoLine.factory(energybounds = ("Fe",2,8)))
        o.markabsorber("Ca",shells=xraylib.K_SHELL,fluolines=element.FluoLine.factory(energybounds = ("Ca",2,8)))
        
        energy = np.array([8,8.5])
        o._fluorescence(energy)
        
        
        
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_fp("test_expi"))
    testSuite.addTest(test_fp("test_transmission"))
    testSuite.addTest(test_fp("test_primary_fluorescence"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

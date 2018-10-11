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
from .. import qxrf
from ...utils import units

import numpy as np

class test_qxrf(unittest.TestCase):

    def geometryinstance(self):
        energy = 10
        geometryinstance = qxrf.factory("QXRFGeometry",diodeI0="iodet1",diodeIt="idet",optics="kb",\
                            xrfdetector="leia",xrfgeometry="sxm120",instrument="sxm",simplecalibration=False)
        
        info = {"I0_counts":300,"It_counts":30,"time":1,"dark":True,"gaindiodeI0":1e8,"gaindiodeIt":1e7}
        geometryinstance.calibrate(**info)
        
        info["I0_counts"] = 10000
        info["It_counts"] = 100000
        info["energy"] = energy-2
        info["dark"] = False
        geometryinstance.calibrate(**info)
        
        info["I0_counts"] = 5000
        info["energy"] = energy+2
        geometryinstance.calibrate(**info)
        
        geometryinstance.setreferenceflux(1e9)
        geometryinstance.setdefaulttime(0.1)
        return geometryinstance
        
    def test_flux(self):
        geometryinstance = self.geometryinstance()

        energy = 10
        time = 0.2
        refflux = 1e9
        
        flux = np.linspace(1e9,1e8,20) # ph/sec
        iodet = geometryinstance.fluxtocps(energy,flux)*time
        
        flux2 = geometryinstance.responsetoflux(energy,iodet/time)
        np.testing.assert_allclose(flux,flux2)
        
        # Normalize data to the real flux (use flux reference)
        rates = np.random.poisson(np.full_like(flux,100)) # 1/ph/sec
        data = flux*time*rates # measured xrf
        dataref = refflux*time*rates # measured when flux whould have been refflux at each poi1e9
        ref = units.Quantity(refflux,"hertz")

        op,_,_,_ = geometryinstance.xrfnormop(energy,expotime=time,reference=ref)
        np.testing.assert_allclose(dataref,data/op(iodet))
        
        # Normalize data to the real flux (use iodet reference)
        iodetref = geometryinstance.fluxtocps(energy,refflux)*time
        ref = units.Quantity(iodetref,"dimensionless")

        op,_,_,_ = geometryinstance.xrfnormop(energy,expotime=time,reference=ref)
        np.testing.assert_allclose(dataref,data/op(iodet))

def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_qxrf("test_flux"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

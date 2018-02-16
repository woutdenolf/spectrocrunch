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
from ..id21_quant import FluxMonitor
from ...common import units

import numpy as np

class test_id21_quant(unittest.TestCase):

    def fluxmonitor(self):
        energy = 10
        monitor = FluxMonitor(iodetname="iodet1",focussed=True,xrfdetector="leia",xrfgeometry="sxm120")
        monitor.setdark(300,None,gainiodet=1e8)
        monitor.setcalib(energy-2,0.5,gainiodet=1e8)
        monitor.setcalib(energy+2,0.5,gainiodet=1e8)
        monitor.setreferenceflux(1e9)
        monitor.settime(0.1)
        return monitor
        
    def test_monitor(self):
        fluxmonitor = self.fluxmonitor()

        energy = 10
        time = 0.2
        refflux = 1e9
        
        flux = np.linspace(1e9,1e8,20) # ph/sec
        iodet = fluxmonitor.fluxtocps(energy,flux)*time
        
        flux2 = fluxmonitor.cpstoflux(energy,iodet/time)
        np.testing.assert_allclose(flux,flux2)
        
        # Normalize data to the real flux (use flux reference)
        rates = np.random.poisson(np.full_like(flux,100)) # 1/ph/sec
        data = flux*time*rates # measured xrf
        dataref = refflux*time*rates # measured when flux whould have been refflux at each poi1e9
        ref = units.Quantity(refflux,"hertz")

        op,_,_,_ = fluxmonitor.xrfnormop(energy,time=time,ref=ref)
        np.testing.assert_allclose(dataref,data/op(iodet))
        
        # Normalize data to the real flux (use iodet reference)
        iodetref = fluxmonitor.fluxtocps(energy,refflux)*time
        ref = units.Quantity(iodetref,"dimensionless")

        op,_,_,_ = fluxmonitor.xrfnormop(energy,time=time,ref=ref)
        np.testing.assert_allclose(dataref,data/op(iodet))

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_id21_quant("test_monitor"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

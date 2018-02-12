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

from .. import xrf as xrfdetectors
from ...geometries import xrf as xrfgeometries
from ...geometries import source as xrfsources

import numpy as np
import silx.math.fit
import scipy.integrate as integrate

class test_xrfdet(unittest.TestCase):                

    def test_lineprofile(self):
        u = 10
        xmin = 0. # not negative, otherwise assert fails (due to step)
        xmax = 2*u
        x = np.linspace(xmin,xmax,1000)
        
        kwargs = {}
        kwargs["mcazero"] = 0. # keV
        kwargs["mcagain"] = 5e-3 # keV
        kwargs["mcanoise"] = 50e-3 # keV
        kwargs["mcafano"] = 0.19
        kwargs["ehole"] = 3.8
        kwargs["activearea"] = 80e-2*0.75 # cm^2
        
        kwargs["shape_conversionenergy"] = u
        
        kwargs["shape_fixedratios"] = {"tailbroadening":0.5,"tailfraction":0.05,"stepfraction":0.005}
        detector1 = xrfdetectors.Detector(**kwargs)
        
        kwargs.pop("shape_fixedratios")
        kwargs["shape_variableratios"] = {"stepheight_ratio":detector1.ratios[1],"tailarea_ratio":detector1.ratios[0],"tailslope_ratio":detector1.tailslope_ratio}
        detector2 = xrfdetectors.Detector(**kwargs)
        
        source = xrfsources.factory("synchrotron")
        
        geometry1 = xrfgeometries.factory("sxm120",detector=detector1,source=source)
        geometry2 = xrfgeometries.factory("sxm120",detector=detector2,source=source)
        
        np.testing.assert_allclose(detector1.fractions,detector2.fractions)
        np.testing.assert_allclose(detector1.ratios,detector2.ratios)
        np.testing.assert_allclose(detector1.tailbroadening,detector2.tailbroadening)
        np.testing.assert_allclose(detector1.tailslope_ratio,detector2.tailslope_ratio)
        
        tmp = str(detector1)
        tmp = str(detector2)
        
        for voigt in [False,True]:
            for tailbroadening in [2.5]: # assert fails when too large
                for wtail in [0,1,0.2,0.4,0.6]:
                    for wstep in [0,1,0.2,0.4,0.6]:
                    

                        # Spectrocrunch arguments
                        if wstep+wtail>1:
                            continue

                        detector1.fractions = (wtail,wstep)
                        detector1.tailbroadening = tailbroadening
                        detector2.ratios = detector1.ratios
                        detector2.tailslope_ratio = detector1.tailslope_ratio

                        if voigt:
                            linewidth = 0.010
                        else:
                            linewidth = 0
                            
                        # Silx arguments
                        wdet = 1-wtail-wstep
                        kwargs = {"gaussian_term":wdet>0,"st_term":wtail>0, "lt_term":False, "step_term":wstep>0}
                        if wdet<=0:
                            garea = 1
                        else:
                            garea = wdet
                        
                        tailarea_ratio,stepheight_ratio = detector1.ratios
                        tailslope_ratio = detector1.tailslope_ratio
                        args = (garea, u, detector1.gaussianFWHM(u), tailarea_ratio, tailslope_ratio, 0, 1, stepheight_ratio)
                        
                        # Compare
                        y1 = np.squeeze(detector1.lineprofile(x,u))
                        y2 = silx.math.fit.sum_ahypermet(x,*args,**kwargs)
                        if wdet>0:
                            y3 = np.squeeze(detector2.lineprofile(x,u))

                        #import matplotlib.pyplot as plt
                        #plt.plot(x,y1+1,'o')
                        #plt.plot(x,y3+1)
                        #plt.show()

                        a = np.max(y1)*0.1
                        np.testing.assert_allclose(y1+a,y2+a)
                        if wdet>0:
                            np.testing.assert_allclose(y1+a,y3+a)
                        
                        # Check unit area
                        area1,error1 = integrate.quad(lambda x: detector1.lineprofile(x,u,linewidth=linewidth),xmin,xmax)
                        area2,error2 = integrate.quad(lambda x: silx.math.fit.sum_ahypermet(np.asarray([x]),*args,**kwargs)[0],xmin,xmax)
                        if voigt:
                            rtol = 1e-3
                        else:
                            rtol = 1e-7
                        np.testing.assert_allclose(area1,1,rtol=rtol)
                        np.testing.assert_allclose(area2,1,rtol=1e-7)
                        if wdet>0:
                            area3,error3 = integrate.quad(lambda x: detector2.lineprofile(x,u,linewidth=linewidth),xmin,xmax)
                            np.testing.assert_allclose(area3,1,rtol=rtol)




def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_xrfdet("test_lineprofile"))
    
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

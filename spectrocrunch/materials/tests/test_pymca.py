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

from .. import multilayer
from .. import pymca
from ...resources import resource_filename
from ...geometries import xrf as xrfgeometries
from ...geometries import source as xraysources
from ...detectors import xrf as xrfdetectors

import numpy as np


class test_pymca(unittest.TestCase):

    def test_rates(self):
        y = np.load(resource_filename("test/mca.npy"))
        cfgfile = resource_filename("test/mca.cfg")
        
        source = xraysources.factory("synchrotron")
        detector = xrfdetectors.factory("XRFDetector",ehole=3.8)
        geometry = xrfgeometries.factory("LinearXRFGeometry",detector=detector,source=source,\
                                    zerodistance=0,detectorposition=0,positionunits="mm")
        sample = multilayer.Multilayer(geometry=geometry)

        pymcahandle = pymca.PymcaHandle(sample=sample,ninteractions=1)
        pymcahandle.loadfrompymca(cfgfile)

        pymcahandle.setdata(y)
        info = pymcahandle.fit()

        # Mass fractions are calculated as follows:
        #   grouparea = flux.time.grouprate
        #   grouprate = solidangle/(4.pi).sum_l[massfrac_l.grouprate_l]  where l loops over the layers
        #
        #   massfrac_l = self.mcafit._fluoRates[layer][element]["rates"][group]      where layer in 1,...,n
        #   grouprate_l = self.mcafit._fluoRates[layer][element]["mass fraction"]
        #   sum_l[massfrac_l.grouprate_l] = self.mcafit._fluoRates[0][element]["rates"][group]*
        #                                   self.mcafit._fluoRates[0][element]["mass fraction"]
        #
        # When element in one layer:
        #   massfrac = area/(flux.time.solidangle/(4.pi).grouprate_l)
        #
        # When element in more than one layer:
        #   grouprate_avg = solidangle/(4.pi).massfrac_avg.sum_l[grouprate_l]
        #
        # When element in more than one layer (per layer as if all intensity came from that layer?):
        #   massfrac_l = grouparea/(flux.time.solidangle/(4.pi).grouprate_l)  

        grouprates = pymcahandle.xraygrouprates(scattering=False,method="fisx")
        safrac = pymcahandle.sample.geometry.solidangle/(4*np.pi)

        for group in info["fitareas"]:
            element,shell = group.split("-")

            grouprate_avg = pymcahandle.mcafit._fluoRates[0][element]["rates"]["{} xrays".format(shell)]
            grouprate_avg *= safrac
            massfrac_avg = pymcahandle.mcafit._fluoRates[0][element]["mass fraction"]

            grouprate = 0.
            grouprate_avg2 = 0.
            massfrac_avg2 = 0.
            npresent = 0
            for j in range(pymcahandle.sample.nlayers):
                i = j+1
                if element in pymcahandle.mcafit._fluoRates[i]:
                    massfrac_l = pymcahandle.mcafit._fluoRates[i][element]["mass fraction"]
                    grouprate_l = pymcahandle.mcafit._fluoRates[i][element]["rates"]["{} xrays".format(shell)]
                    grouprate += massfrac_l*grouprate_l
                    grouprate_avg2 += grouprate_l
                    massfrac_avg2 = max(massfrac_avg2,massfrac_l)
                    npresent += 1
            grouprate *= safrac
            grouprate_avg2 *= safrac
            if group in grouprates:
                np.testing.assert_allclose(grouprate,grouprates[group],rtol=1e-2) # 1% error fisx vs. Elements?
            
            np.testing.assert_allclose(grouprate,info["rates"][group])
            
            if npresent==1:
                np.testing.assert_allclose(grouprate_avg*massfrac_avg,grouprate)
            else:
                np.testing.assert_allclose(massfrac_avg,massfrac_avg2)
                np.testing.assert_allclose(grouprate_avg,grouprate_avg2)

            np.testing.assert_allclose(info["massfractions"][group],info["fitareas"][group]/(pymcahandle.flux*pymcahandle.time*grouprate_avg))

            for j in range(pymcahandle.sample.nlayers):
                i = j+1
                if element in pymcahandle.mcafit._fluoRates[i]:
                    grouprate_l = pymcahandle.mcafit._fluoRates[i][element]["rates"]["{} xrays".format(shell)]
                    grouprate = grouprate_l*safrac
                    np.testing.assert_allclose(info["lmassfractions"][j][group],info["fitareas"][group]/(pymcahandle.flux*pymcahandle.time*grouprate))

        if False:
            import matplotlib.pyplot as plt
            plt.plot(data["x"],data["y"],label='data')
            plt.plot(data["x"],data["yfit"],label='pymca')
            
            spectrum = pymcahandle.xrayspectrum()
            spectrum.plot(fluxtime=pymcahandle.flux*pymcahandle.time,histogram=True,log=False,decompose=False)

            ax = plt.gca()
            ax.set_ylim(ymin=np.nanmin(y[np.nonzero(y)]))
            ax.set_xlabel("Energy (keV)")
            ax.set_ylabel("Intensity (cts)")
            plt.legend(loc='best')
            plt.show()
        
        
def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_pymca("test_rates"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
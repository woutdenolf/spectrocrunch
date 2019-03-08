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
import numpy as np
import silx.math.fit
import scipy.integrate as integrate
import itertools

from .. import xrf
from ...patch import jsonpickle


class test_xrf(unittest.TestCase):

    def plot(self, x, y1, y2):
        import matplotlib.pyplot as plt
        plt.plot(x, y1+1)
        plt.plot(x, y2+1)
        plt.show()

    def test_lineprofile(self):
        u = 10
        xmin = 0.  # not negative, otherwise assert fails (due to step)
        xmax = 2*u
        x = np.linspace(xmin, xmax, 1000)

        kwargs = {}
        kwargs["mcazero"] = 0.  # keV
        kwargs["mcagain"] = 5e-3  # keV
        kwargs["mcanoise"] = 50e-3  # keV
        kwargs["mcafano"] = 0.19
        kwargs["ehole"] = 3.8
        kwargs["activearea"] = 80e-2*0.75  # cm^2

        kwargs["shape_conversionenergy"] = u

        kwargs["shape_fixedarearatios"] = {"stailbroadening": 0.1, "stailfraction": 0.05,
                                           "ltailbroadening": 0.5, "ltailfraction": 0.01,
                                           "stepfraction": 0.005,
                                           "bpeak": True, "bstail": True, "bltail": True, "bstep": True}
        detector1 = xrf.XRFDetector(**kwargs)

        kwargs.pop("shape_fixedarearatios")
        kwargs["shape_pymca"] = {"stailarea_ratio": detector1.ratios[0],
                                 "stailslope_ratio": detector1.stailslope_ratio,
                                 "ltailarea_ratio": detector1.ratios[1],
                                 "ltailslope_ratio": detector1.ltailslope_ratio,
                                 "stepheight_ratio": detector1.ratios[2],
                                 "bpeak": True, "bstail": True, "bltail": True, "bstep": True}
        detector2 = xrf.XRFDetector(**kwargs)

        kwargs.pop("shape_pymca")
        kwargs["shape_fixedarearatios"] = {"stailbroadening": detector1.stailbroadening,
                                           "stailfraction": detector1.fractions[0],
                                           "ltailbroadening": detector1.ltailbroadening,
                                           "ltailfraction": detector1.fractions[1],
                                           "stepfraction": detector1.fractions[2],
                                           "bpeak": True, "bstail": True, "bltail": True, "bstep": True}
        detector3 = xrf.XRFDetector(**kwargs)

        for a, b in itertools.permutations([detector1, detector2, detector3], 2):
            np.testing.assert_allclose(a.all_fractions, b.all_fractions)
            np.testing.assert_allclose(a.ratios, b.ratios)
            np.testing.assert_allclose(a.stailbroadening, b.stailbroadening)
            np.testing.assert_allclose(a.stailslope_ratio, b.stailslope_ratio)
            np.testing.assert_allclose(a.ltailbroadening, b.ltailbroadening)
            np.testing.assert_allclose(a.ltailslope_ratio, b.ltailslope_ratio)

        tmp = str(detector1)
        tmp = str(detector2)
        tmp = str(detector2)

        parameters = [[False, True],
                      [2.5],  # assert when too large
                      [1.5],  # assert when too large
                      [0, 1, 0.2, 0.4, 0.6],
                      [0, 1, 0.2, 0.4, 0.6],
                      [0, 1, 0.2, 0.4, 0.6],
                      [False, True]]
        for params in itertools.product(*parameters):
            voigt, ltailbroadening, stailbroadening, wstep, wstail, wltail, normalized = params
            # Spectrocrunch arguments
            try:
                wpeak = detector1.wpeak(wstail, wltail, wstep)
            except RuntimeError:
                continue
            bpeak = wpeak > 0
            bstail = wstail > 0
            bltail = wltail > 0
            bstep = wstep > 0
            detector1.bpeak = bpeak
            detector2.bpeak = bpeak
            detector1.bstail = bstail
            detector1.bltail = bltail
            detector2.bstail = bstail
            detector2.bltail = bltail
            detector1.bstep = bstep
            detector2.bstep = bstep

            detector1.stailbroadening = stailbroadening
            detector1.ltailbroadening = ltailbroadening
            detector1.fractions = (wstail, wltail, wstep)

            detector2.stailslope_ratio = detector1.stailslope_ratio
            detector2.ltailslope_ratio = detector1.ltailslope_ratio
            detector2.ratios = detector1.ratios

            # Silx arguments
            kwargs = {"gaussian_term": wpeak > 0, "st_term": wstail >
                      0, "lt_term": wltail > 0, "step_term": wstep > 0}
            if bpeak and normalized:
                garea = wpeak
            else:
                garea = 1
            stailarea_ratio, ltailarea_ratio, stepheight_ratio = detector1.ratios
            stailslope_ratio = detector1.stailslope_ratio
            ltailslope_ratio = detector1.ltailslope_ratio
            args = (garea, u, detector1.gaussianFWHM(u),
                    stailarea_ratio, stailslope_ratio,
                    ltailarea_ratio, ltailslope_ratio,
                    stepheight_ratio)

            # Compare
            y1 = np.squeeze(detector1.lineprofile(x, u, normalized=normalized))
            y2 = silx.math.fit.sum_ahypermet(x, *args, **kwargs)
            y3 = np.squeeze(detector2.lineprofile(x, u, normalized=normalized))

            a = np.max(y1)*0.1
            # print wpeak,wstail,wltail,wstep
            # self.plot(x,y2,y3)
            rtol = 1e-6
            np.testing.assert_allclose(y1+a, y2+a, rtol=rtol)
            np.testing.assert_allclose(y2+a, y3+a, rtol=rtol)

            # Check unit area
            if voigt:
                linewidth = 0.010
                rtol = 1e-3
            else:
                linewidth = 0
                rtol = 1e-7
            area1, error1 = integrate.quad(lambda x: detector1.lineprofile(
                x, u, linewidth=linewidth, normalized=normalized), xmin, xmax)
            area2, error2 = integrate.quad(lambda x: silx.math.fit.sum_ahypermet(
                np.asarray([x]), *args, **kwargs)[0], xmin, xmax)
            area3, error3 = integrate.quad(lambda x: detector2.lineprofile(
                x, u, linewidth=linewidth, normalized=normalized), xmin, xmax)
            if normalized:
                np.testing.assert_allclose(area1, 1, rtol=rtol)
                np.testing.assert_allclose(area2, 1, rtol=1e-7)
                np.testing.assert_allclose(area3, 1, rtol=rtol)
            # else:
            #    print area1,area2,area3
            #    np.testing.assert_raises(AssertionError, np.testing.assert_allclose, area1,1,rtol=rtol)
            #    np.testing.assert_raises(AssertionError, np.testing.assert_allclose, area2,1,rtol=1e-7)
            #    np.testing.assert_raises(AssertionError, np.testing.assert_allclose, area3,1,rtol=rtol)

    @unittest.skipIf(xrf.compoundfromname.xraylib is None,
                     "xraylib not installed")
    def test_convert(self):
        exclude = 'XRFDetector',
        for name, cls in xrf.XRFDetector.clsregistry.items():
            if name not in exclude:
                a = cls(shape_conversionenergy=7.)
                b = a.convert(inplace=False)
                b.convert(inplace=True)
                c = b.convert(inplace=False)
                c.convert(inplace=True)
                np.testing.assert_allclose(a.all_fractions, b.all_fractions)
                np.testing.assert_allclose(c.all_fractions, b.all_fractions)
                np.testing.assert_allclose(a.ratios, b.ratios)
                np.testing.assert_allclose(c.ratios, b.ratios)
                np.testing.assert_allclose(a.stailbroadening, b.stailbroadening)
                np.testing.assert_allclose(c.stailbroadening, b.stailbroadening)
                np.testing.assert_allclose(a.ltailbroadening, b.ltailbroadening)
                np.testing.assert_allclose(c.ltailbroadening, b.ltailbroadening)
                np.testing.assert_allclose(a.stailslope_ratio, b.stailslope_ratio)
                np.testing.assert_allclose(c.stailslope_ratio, b.stailslope_ratio)
                np.testing.assert_allclose(a.ltailslope_ratio, b.ltailslope_ratio)
                np.testing.assert_allclose(c.ltailslope_ratio, b.ltailslope_ratio)

    @unittest.skipIf(xrf.compoundfromname.xraylib is None,
                     "xraylib not installed")
    def test_serialize(self):
        exclude = 'XRFDetector',
        for name, cls in xrf.XRFDetector.clsregistry.items():
            if name not in exclude:
                d1 = cls()
                d2 = jsonpickle.loads(jsonpickle.dumps(d1))
                self.assertEqual(d1, d2)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_xrf("test_serialize"))
    testSuite.addTest(test_xrf("test_convert"))
    testSuite.addTest(test_xrf("test_lineprofile"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

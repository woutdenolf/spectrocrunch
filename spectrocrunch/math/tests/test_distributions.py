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

from .. import distributions

from scipy.integrate import quad
import numpy as np
import scipy.stats


class test_distributions(unittest.TestCase):

    def _check_pdfintegral(self, integral, integrale, theory):
        integrale = max(integrale, 1e-5)
        limits = integral+integrale*np.array([-1, 1])
        self.assertTrue(limits[0] <= theory <= limits[1])

    def quad(self, func, a, b):
        return quad(func, a, b)
        if b >= a:
            return quad(func, a, b)
        else:
            return 0., 1e-16

    def func(self, func, args=None):
        if args:
            return lambda x: func(x, *args)
        else:
            return func

    def checkpdf(self, rv, qmin, qmax, xmin, xmax, xn, args=None):
        # self.plot(rv,qmin,qmax,xmin,xmax,xn,args=args)
        # return

        pdf = self.func(rv.pdf, args=args)
        cdf = self.func(rv.cdf, args=args)
        ppf = self.func(rv.ppf, args=args)

        integral, integrale = self.quad(pdf, qmin, qmax)
        self._check_pdfintegral(integral, integrale, 1)

        x = np.linspace(xmin, xmax, xn)
        np.allclose(ppf(cdf(x)), x)

        p = np.linspace(0, 1, xn)
        # TODO: check RuntimeWarning
        np.allclose(cdf(ppf(x)), x)

        for x in np.linspace(xmin, xmax, xn):
            integral, integrale = self.quad(pdf, qmin, x)
            self._check_pdfintegral(integral, integrale, cdf(x))

    def plot(self, rv, qmin, qmax, xmin, xmax, xn, args=None):
        import matplotlib.pyplot as plt

        pdf = self.func(rv.pdf, args=args)
        cdf = self.func(rv.cdf, args=args)
        ppf = self.func(rv.ppf, args=args)

        x = np.linspace(xmin, xmax, xn)
        plt.plot(x, pdf(x), 'o-')
        plt.plot(x, cdf(x), 'o-')
        plt.plot(x, [self.quad(pdf, qmin, xi)[0] for xi in x], 'o-')
        plt.show()

    def test_pdf(self):
        qmin, qmax = -np.inf, np.inf
        xmin, xmax, xn = -10, 10, 50
        rv = scipy.stats.norm
        self.checkpdf(rv, qmin, qmax, xmin, xmax, xn)

        #k = 1
        #qmin,qmax = -k,k
        #xmin,xmax,xn = -10,10,50
        #rv = scipy.stats.truncnorm(a=qmin,b=qmax)
        # self.checkpdf(rv,qmin,qmax,xmin,xmax,xn)

        k = 1
        qmin, qmax = -k, k
        xmin, xmax, xn = -10, 10, 50
        rv = distributions.limitednorm(k)
        self.checkpdf(rv, qmin, qmax, xmin, xmax, xn)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_distributions("test_pdf"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

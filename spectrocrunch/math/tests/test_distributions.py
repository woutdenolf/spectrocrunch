# -*- coding: utf-8 -*-

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

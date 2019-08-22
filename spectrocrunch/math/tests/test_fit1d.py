# -*- coding: utf-8 -*-

import unittest

from .. import fit1d

import numpy as np


class test_fit1d(unittest.TestCase):

    def test_leastsq(self):
        nx = 501
        x = np.arange(nx)
        x0 = 10
        sx = nx//4
        A = 1000.
        p1 = np.array([x0, sx, A], dtype=np.float32)
        x0, sx, A = tuple(p1)

        data = fit1d.gaussian(x, x0, sx, A)
        #import matplotlib.pyplot as plt
        # plt.plot(data)
        # plt.show()

        p2, _ = fit1d.fitgaussian(x, data)
        np.testing.assert_allclose(p1, p2)

    def test_linfit(self):
        x = np.asarray([1.47, 1.50, 1.52, 1.55, 1.57, 1.60, 1.63,
                        1.65, 1.68, 1.70, 1.73, 1.75, 1.78, 1.80, 1.83])
        y = np.asarray([52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 59.93,
                        61.29, 63.11, 64.47, 66.28, 68.10, 69.92, 72.19, 74.46])
        m = 61.272
        b = -39.062
        vm = 3.1539
        vb = 8.63185
        (m1, b1), (em1, eb1) = fit1d.linfit(x, y, errors=True)
        (m2, b2), (em2, eb2) = fit1d.linfit2(x, y, errors=True)
        vm1 = em1*em1
        vm2 = em2*em2
        vb1 = eb1*eb1
        vb2 = eb2*eb2
        np.testing.assert_array_almost_equal([m, b], [m1, b1], decimal=3)
        np.testing.assert_array_almost_equal([m, b], [m2, b2], decimal=3)
        np.testing.assert_array_almost_equal([vm, vm], [vm1, vm2], decimal=3)
        np.testing.assert_array_almost_equal([vb, vb], [vb1, vb2], decimal=4)
        np.testing.assert_allclose([m1, b1], [m2, b2])

        m1, em1 = fit1d.linfit_zerointercept(x, y, errors=True)
        m2, em2 = fit1d.linfit_zerointercept2(x, y, errors=True)
        np.testing.assert_allclose(m1, m2)
        np.testing.assert_almost_equal(
            em1, em2, decimal=3)  # em2 maybe not correct???


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_fit1d("test_leastsq"))
    testSuite.addTest(test_fit1d("test_linfit"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

# -*- coding: utf-8 -*-

import unittest

from .. import fit2d

import numpy as np
import pylab


class test_fit2d(unittest.TestCase):
    def test_leastsq(self):
        nx = 501
        ny = 401
        y, x = np.indices((ny, nx))
        x0 = 10
        y0 = ny // 2
        sx = nx // 4
        sy = ny // 4
        rho = 0.5
        A = 1000.0
        p1 = np.array([x0, y0, sx, sy, rho, A], dtype=np.float32)
        x0, y0, sx, sy, rho, A = tuple(p1)

        data = fit2d.gaussian(x, y, x0, y0, sx, sy, rho, A)
        # self.plot(data)

        p2, _ = fit2d.fitgaussian(x, y, data)
        np.testing.assert_allclose(p1, p2)

    def plot(self, img):
        pylab.figure(1)
        pylab.subplot(111)
        pylab.imshow(img, origin="lower", interpolation="nearest")
        pylab.pause(0.1)
        raw_input("Press enter to continue...")


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_fit2d("test_leastsq"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

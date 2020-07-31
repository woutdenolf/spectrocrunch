# -*- coding: utf-8 -*-

import unittest
import numpy as np
import itertools

from .. import interpolate


class test_interpolate(unittest.TestCase):
    def _generate(self, shape):
        axes = [range(n) for n in shape]
        data = np.random.normal(size=shape).astype(np.float32) * 100
        return data, axes

    def test_nd(self):
        def identity(x):
            return x

        def meshgrid(axes):
            return np.meshgrid(*axes, indexing="ij")

        parameters = (
            (0, 1, 2),
            (True,),
            (
                (interpolate.interpolate_regular, identity),
                (interpolate.interpolate_irregular, meshgrid),
            ),
        )

        for parameters in itertools.product(*parameters):
            degree, asgrid, (interp, meshgrid) = parameters

            rtol = 1e-3

            # 1D
            data1, axes = self._generate((6,))
            axnew = (0,)
            data2 = interp(data1, axes, axnew, degree=degree, asgrid=asgrid)
            data3 = data1[axnew]
            np.testing.assert_allclose(data2, data3, rtol=rtol)
            axnew = ([0, 5],)
            data2 = interp(data1, axes, axnew, degree=degree, asgrid=asgrid)
            data3 = data1[axnew]
            np.testing.assert_allclose(data2, data3, rtol=rtol)

            # 2D
            data1, axes = self._generate((6, 8))
            if asgrid:
                axnew = ([0, 3, 5], [1, 3, 6])
            else:
                axnew = ([0, 3, 5], [1, 3])
            data2 = interp(data1, meshgrid(axes), axnew, degree=degree, asgrid=asgrid)
            if asgrid:
                axnew = tuple(np.meshgrid(*axnew, indexing="ij"))
            data3 = data1[axnew]
            np.testing.assert_allclose(data2, data3, rtol=rtol)

            # 3D
            data1, axes = self._generate((6, 8, 9))
            if asgrid:
                axnew = ([0, 3, 5], [1, 3, 6], [2, 5, 7])
            else:
                axnew = ([0, 3, 5], [1, 3], [2, 7])
            data2 = interp(data1, meshgrid(axes), axnew, degree=degree, asgrid=asgrid)
            if asgrid:
                axnew = tuple(np.meshgrid(*axnew, indexing="ij"))
            data3 = data1[axnew]
            np.testing.assert_allclose(data2, data3, rtol=rtol)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_interpolate("test_nd"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

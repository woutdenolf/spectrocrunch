# -*- coding: utf-8 -*-

import unittest
import numpy as np
from scipy.stats import ortho_group

from .. import quadrics


class test_quadrics(unittest.TestCase):
    def _assert_matrix(self, eval_matrix, eval_eqn, nargs, npts=100):
        args = [np.random.uniform(-10, 10, (npts, 3)) for _ in range(nargs)]
        xT = np.random.uniform(-10, 10, (npts, 3))
        ones = np.ones((npts, 1))
        xhT = np.hstack([xT, ones])
        for argsi in zip(*args):
            A = eval_matrix(*argsi)
            actual = (xhT.dot(A) * xhT).sum(axis=-1)
            argsi = [np.atleast_2d(a) for a in argsi]
            expected = eval_eqn(xT, *argsi)
            np.testing.assert_allclose(actual, expected)

    def test_plane(self):
        def eval_matrix(x0, u):
            return quadrics.plane(x0, u)

        def eval_eqn(xT, x0, u):
            return ((xT - x0) * u).sum(axis=-1)

        self._assert_matrix(eval_matrix, eval_eqn, 2)

    def test_cylinder(self):
        def eval_matrix(x0, u):
            x0 = x0[[0, 2]]
            u = u[[0, 2]]
            return quadrics.cylinder(x0, u, 1)

        def eval_eqn(xT, x0, u):
            xT = xT[:, [0, 2]]
            x0 = x0[:, [0, 2]]
            u = u[:, [0, 2]]
            return ((xT - x0) ** 2 / u ** 2).sum(axis=1) - 1

        self._assert_matrix(eval_matrix, eval_eqn, 2)

    def test_ellipsoid(self):
        def eval_matrix(x0, u):
            return quadrics.ellipsoid(x0, u)

        def eval_eqn(xT, x0, u):
            return ((xT - x0) ** 2 / u ** 2).sum(axis=1) - 1

        self._assert_matrix(eval_matrix, eval_eqn, 2)

    def test_transform_plane(self):
        def to_cart(x):
            return x[:-1] / x[-1]

        def to_hom(x):
            return np.append(x, 1)

        C = np.eye(4)
        x0, u, C[0:3, 3] = np.random.uniform(-10, 10, (3, 3))
        C[0:3, 0:3] = ortho_group.rvs(3)
        L = np.linalg.inv(C)

        x0_2 = to_cart(L.dot(to_hom(x0)))
        u_2 = L[0:3, 0:3].dot(x0 + u) - L[0:3, 0:3].dot(x0)

        A = quadrics.plane(x0, u)
        A_2 = quadrics.plane(x0_2, u_2)
        A_3 = quadrics.change_of_frame(A, C)

        np.testing.assert_allclose(A_2, A_3)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_quadrics("test_plane"))
    testSuite.addTest(test_quadrics("test_cylinder"))
    testSuite.addTest(test_quadrics("test_ellipsoid"))
    testSuite.addTest(test_quadrics("test_transform_plane"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

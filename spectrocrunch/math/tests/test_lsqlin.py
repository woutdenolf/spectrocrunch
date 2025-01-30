# -*- coding: utf-8 -*-

import unittest
import numpy as np

from .. import lsqlin


class test_lsqlin(unittest.TestCase):
    def test_lsqlin(self):

        C = np.array(
            [
                [0.9501, 0.7620, 0.6153, 0.4057],
                [0.2311, 0.4564, 0.7919, 0.9354],
                [0.6068, 0.0185, 0.9218, 0.9169],
                [0.4859, 0.8214, 0.7382, 0.4102],
                [0.8912, 0.4447, 0.1762, 0.8936],
            ]
        )
        sC = lsqlin.sparse.coo_matrix(C)
        csC = lsqlin.scipy_sparse_to_spmatrix(sC)

        A = np.array(
            [
                [0.2027, 0.2721, 0.7467, 0.4659],
                [0.1987, 0.1988, 0.4450, 0.4186],
                [0.6037, 0.0152, 0.9318, 0.8462],
            ]
        )
        sA = lsqlin.sparse.coo_matrix(A)
        csA = lsqlin.scipy_sparse_to_spmatrix(sA)

        d = np.array([0.0578, 0.3528, 0.8131, 0.0098, 0.1388])
        md = lsqlin.matrix(d)

        b = np.array([0.5251, 0.2026, 0.6721])
        mb = lsqlin.matrix(b)

        lb = np.array([-0.1] * 4)
        mlb = lsqlin.matrix(lb)
        mmlb = -0.1

        ub = np.array([2] * 4)
        mub = lsqlin.matrix(ub)
        mmub = 2

        opts = {"show_progress": False}
        sol = [-1.00e-01, -1.00e-01, 2.15e-01, 3.50e-01]

        for iC in [C, sC, csC]:
            for iA in [A, sA, csA]:
                for iD in [d, md]:
                    for ilb in [lb, mlb, mmlb]:
                        for iub in [ub, mub, mmub]:
                            for ib in [b, mb]:
                                ret = lsqlin.lsqlin(
                                    iC, iD, reg=0, A=iA, b=ib, lb=ilb, ub=iub, opts=opts
                                )
                                result = lsqlin.cvxopt_to_numpy_matrix(ret["x"])
                                np.testing.assert_allclose(result, sol, rtol=2e-3)

    def test_lsqnonneg(self):
        C = np.array(
            [[0.0372, 0.2869], [0.6861, 0.7071], [0.6233, 0.6245], [0.6344, 0.6170]]
        )
        d = np.array([0.8587, 0.1781, 0.0747, 0.8405])
        ret = lsqlin.lsqnonneg(C, d)
        result = lsqlin.cvxopt_to_numpy_matrix(ret["x"])
        sol = [2.5e-07, 6.93e-01]
        np.testing.assert_allclose(result, sol, rtol=2e-3)


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_lsqlin("test_lsqlin"))
    testSuite.addTest(test_lsqlin("test_lsqnonneg"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

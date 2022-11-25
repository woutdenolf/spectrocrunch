# -*- coding: utf-8 -*-

import unittest
import itertools

from .. import noisepropagation
from .. import fit1d
from ...utils import units
from ...utils import instance

import numpy as np
import uncertainties


class test_noisepropagation(unittest.TestCase):
    def test_instance(self):
        for RV, Q, expand, arr in itertools.product([True, False], repeat=4):
            if arr:
                x = np.arange(1.0, 10)
            else:
                x = 10.0

            if RV:
                s = np.sqrt(x)
                y = noisepropagation.randomvariable(x, s)
            else:
                s = x * 0
                y = x

            if Q:
                unit = "mm"
                y = units.Quantity(y, units=unit)
            else:
                unit = None
            if expand:
                y = instance.asarray(y)

            np.testing.assert_array_equal(
                x, units.magnitude(noisepropagation.E(y), units=unit)
            )
            np.testing.assert_array_equal(s, noisepropagation.S(y))

        x = 10, noisepropagation.randomvariable(20, 1)
        np.testing.assert_array_equal([10, 20], noisepropagation.E(x))
        np.testing.assert_array_equal([0, 1], noisepropagation.S(x))

    def test_repeat(self):
        X = noisepropagation.randomvariable([100, 200], [100**0.5, 200**0.5])

        N = 10
        Y1 = noisepropagation.repeat(N, X)
        Y2 = noisepropagation.compound(noisepropagation.randomvariable(N, 0), X)
        Y2 = np.squeeze(Y2)

        np.testing.assert_allclose(noisepropagation.E(Y1), noisepropagation.E(Y2))
        np.testing.assert_allclose(noisepropagation.S(Y1), noisepropagation.S(Y2))
        np.testing.assert_allclose(noisepropagation.SNR(X), noisepropagation.SNR(N * X))
        np.testing.assert_allclose(
            noisepropagation.SNR(Y1), noisepropagation.SNR(X) * N**0.5
        )

    def _RVequal(self, N1, N2):
        self.assertEqual(noisepropagation.E(N1), noisepropagation.E(N2))
        self.assertEqual(noisepropagation.S(N1), noisepropagation.S(N2))

    def _RValmostequal(self, N1, N2):
        np.testing.assert_allclose(noisepropagation.E(N1), noisepropagation.E(N2))
        np.testing.assert_allclose(noisepropagation.S(N1), noisepropagation.S(N2))

    def _RVnotequal(self, N1, N2):
        self.assertNotEqual(noisepropagation.E(N1), noisepropagation.E(N2))
        self.assertNotEqual(noisepropagation.S(N1), noisepropagation.S(N2))

    def test_compound(self):
        N0 = noisepropagation.randomvariable([100.0, 200], [100**0.5, 200**0.5])
        EN0 = noisepropagation.E(N0)

        g1 = noisepropagation.randomvariable(10.0, 1.3)
        g2 = noisepropagation.randomvariable(3.0, 0.2)
        g3 = noisepropagation.randomvariable(7.0, 1.5)
        g4 = noisepropagation.randomvariable(6.0, 1.7)

        Eg1 = noisepropagation.E(g1)
        Eg2 = noisepropagation.E(g2)
        Eg3 = noisepropagation.E(g3)
        Eg4 = noisepropagation.E(g4)

        N1 = noisepropagation.compound(N0, g1)
        N2 = noisepropagation.compound(N1, g2)
        N3 = noisepropagation.compound(N2, g3)
        N4 = noisepropagation.compound(N3, g4)

        a = Eg1 * Eg2 * Eg3 * Eg4

        b = (
            noisepropagation.VAR(g1) * (Eg2 * Eg3 * Eg4) ** 2
            + noisepropagation.VAR(g2) * (Eg1) * (Eg3 * Eg4) ** 2
            + noisepropagation.VAR(g3) * (Eg1 * Eg2) * (Eg4) ** 2
            + noisepropagation.VAR(g4) * (Eg1 * Eg2 * Eg3)
        )
        np.testing.assert_allclose(noisepropagation.E(N4), a * EN0)
        np.testing.assert_allclose(
            noisepropagation.VAR(N4), b * EN0 + a**2 * noisepropagation.VAR(N0)
        )

        np.testing.assert_allclose(
            noisepropagation.RVAR(N4),
            noisepropagation.RVAR(N0)
            + (
                noisepropagation.RVAR(g1)
                + noisepropagation.RVAR(g2) / Eg1
                + noisepropagation.RVAR(g3) / (Eg1 * Eg2)
                + noisepropagation.RVAR(g4) / (Eg1 * Eg2 * Eg3)
            )
            / EN0,
        )

    def test_bernouilli(self):
        # Bernouilli processes can be multiplied
        g1 = noisepropagation.bernouilli(0.5)
        g2 = noisepropagation.bernouilli(0.4)
        g3 = noisepropagation.bernouilli(
            noisepropagation.E(g1) * noisepropagation.E(g2)
        )

        N = noisepropagation.randomvariable(100.0, 12.3)

        N1 = noisepropagation.compound(noisepropagation.compound(N, g1), g2)
        N2 = noisepropagation.compound(N, g3)
        self._RVequal(N1, N2)

        # Bernouilli compounding of Poisson gives Poisson
        N = noisepropagation.poisson(100)
        N1 = noisepropagation.compound(N, g1)
        N2 = noisepropagation.poisson(noisepropagation.E(N) * noisepropagation.E(g1))
        self._RVequal(N1, N2)

    def test_poisson(self):
        # Poisson processes cannot be multiplied
        g1 = noisepropagation.poisson(10)
        g2 = noisepropagation.poisson(15)
        g3 = noisepropagation.poisson(noisepropagation.E(g1) * noisepropagation.E(g2))

        N = noisepropagation.randomvariable(100.0, 12.3)

        N1 = noisepropagation.compound(noisepropagation.compound(N, g1), g2)
        N2 = noisepropagation.compound(N, g3)
        self.assertEqual(noisepropagation.E(N1), noisepropagation.E(N2))
        self.assertNotEqual(noisepropagation.S(N1), noisepropagation.S(N2))

        # Poisson compounding of Poisson does not gives Poisson
        N = noisepropagation.poisson(100)
        N1 = noisepropagation.compound(N, g1)
        N2 = noisepropagation.poisson(noisepropagation.E(N) * noisepropagation.E(g1))
        self.assertEqual(noisepropagation.E(N1), noisepropagation.E(N2))
        self.assertNotEqual(noisepropagation.S(N1), noisepropagation.S(N2))

    def test_reverse(self):

        X = noisepropagation.randomvariable(10.0, 1.3)
        X2 = noisepropagation.randomvariable(10.0, 1.3)
        Y = noisepropagation.randomvariable(3.0, 0.2)
        Y2 = noisepropagation.randomvariable(3.0, 0.2)

        Z = X + Y
        XX = noisepropagation.reverse_add(Z, Y2)
        self._RValmostequal(X, XX)

        Z = X - Y
        XX = noisepropagation.reverse_sub(Z, Y2)
        self._RValmostequal(X, XX)

        Z = X * Y
        XX = noisepropagation.reverse_mult(Z, Y2)
        self._RValmostequal(X, XX)

        Z = X / Y
        XX = noisepropagation.reverse_div(Z, Y2)
        self._RValmostequal(X, XX)

        Z = uncertainties.unumpy.log(X2)
        XX = noisepropagation.reverse_log(Z)
        self._RValmostequal(X, XX)

        Z = noisepropagation.repeat(10, X)
        XX = noisepropagation.repeat(10, Z, forward=False)
        self._RValmostequal(X, XX)

        Z = noisepropagation.compound(X, Y)
        XX = noisepropagation.compound(Z, Y, forward=False)
        self._RValmostequal(X, XX)

    def test_lstsq_std(self):
        nx, npa = 10, 3
        A = np.random.rand(nx, npa) * 10

        # x -LINPROP-> VAR(b)
        # VAR(x) -LINPROP-> VAR(b)
        x = np.random.rand(npa)
        stdx = np.random.rand(npa)
        x = noisepropagation.randomvariable(x, stdx)
        varb = noisepropagation.VAR(np.dot(A, x))
        varb2 = np.dot(A * A, stdx**2)
        np.testing.assert_allclose(varb, varb2)

        stdx2 = np.sqrt(fit1d.lstsq(A * A, varb))
        np.testing.assert_allclose(stdx, stdx2)

        stdx2 = fit1d.lstsq_std_indep(A, vare=varb)
        np.testing.assert_allclose(stdx, stdx2)
        # stdx3 = fit1d.lstsq_std(A,vare=varb)
        # covx = fit1d.lstsq_cov(A,vare=varb)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_noisepropagation("test_instance"))
    testSuite.addTest(test_noisepropagation("test_repeat"))
    testSuite.addTest(test_noisepropagation("test_compound"))
    testSuite.addTest(test_noisepropagation("test_bernouilli"))
    testSuite.addTest(test_noisepropagation("test_poisson"))
    testSuite.addTest(test_noisepropagation("test_reverse"))
    testSuite.addTest(test_noisepropagation("test_lstsq_std"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

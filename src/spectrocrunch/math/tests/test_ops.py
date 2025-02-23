import unittest
import itertools
import numpy as np
import random

from .. import linop


class test_ops(unittest.TestCase):
    def test_linop(self):
        for x in range(-2, 2):
            for m1 in range(-11, 11, 3):
                for b1 in range(-11, 11, 3):
                    o1 = linop.LinearOperator(m1, b1)
                    o1i = o1.inverse

                    self.assertAlmostEqual(o1(x), o1.m * x + o1.b)
                    self.assertAlmostEqual(o1i(x), (x - o1.b) / float(o1.m))
                    self.assertEqual(o1 * o1, o1**2)
                    self.assertEqual(o1 * o1 * o1, o1**3)
                    self.assertAlmostEqual((o1i * o1i)(x), (o1 * o1).inverse(x))
                    self.assertAlmostEqual(
                        (o1i * o1i * o1i)(x), (o1 * o1 * o1).inverse(x)
                    )

                    for m2 in range(-11, 11, 3):
                        for b2 in range(-11, 11, 3):
                            o2 = linop.LinearOperator(m2, b2)
                            o2i = o2.inverse
                            self.assertAlmostEqual((o1 * o2)(x), o1(o2(x)))
                            self.assertAlmostEqual((o1 * o2).inverse(x), o2i(o1i(x)))

    def _gencase(self, ncases=100):
        sops = [
            linop.LinearOperator(1.3, 0.1),
            linop.Identity(),
            linop.LinearOperator(-0.9, -0.3),
            linop.Clip(-10, 10),
            linop.Clip(None, 10),
            linop.Clip(10, None),
            None,
            linop.NaNClip(-5, 11),
            linop.Clip(3, 7),
            linop.NaNClip(-7, 0),
            linop.NaNClip(None, 0),
            linop.NaNClip(-7, None),
        ]

        for i in range(ncases):
            n = random.randint(1, 5)
            combinations = list(itertools.combinations_with_replacement(sops, n))
            for combination in random.sample(combinations, n):
                yield combination

    def test_combine(self):
        for arg in np.random.random(5) * 60 - 30:
            for ops in self._gencase(ncases=500):
                opc = linop.Identity()
                result = arg
                # Forward
                for op in ops:
                    if op is not None:
                        result = op(result)
                    opc = op * opc
                    np.testing.assert_allclose(result, opc(arg))
                # Inverse
                result = arg
                for op in reversed(ops):
                    if op is not None:
                        result = op.inverse(result)
                np.testing.assert_allclose(result, opc.inverse(arg))

    def test_veclinop(self):
        m = np.linspace(1, 5, 5)[:, np.newaxis]
        b = np.linspace(1, 5, 5)[:, np.newaxis]
        op = linop.LinearOperator(m, b)
        x = np.arange(10)[np.newaxis, :]
        y1 = op(x)
        y2 = m * x + b
        self.assertEqual(y1.shape, (5, 10))
        np.testing.assert_array_equal(y1, y2)
        self.assertEqual(op, op)
        op *= op
        self.assertEqual(op, op)


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_ops("test_linop"))
    testSuite.addTest(test_ops("test_combine"))
    testSuite.addTest(test_ops("test_veclinop"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

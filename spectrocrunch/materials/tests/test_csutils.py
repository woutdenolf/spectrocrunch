# -*- coding: utf-8 -*-

import unittest

from ..import csutils
import numpy as np


class test_csutils(unittest.TestCase):

    @unittest.skipIf(csutils.xraylib.XRayInit is None,
                     "xraylib not installed")
    def test_shape(self):
        keep = csutils.xraylib.xraylib_np
        csutils.xraylib.xraylib_np = None
        try:
            self.assert_shape()
        finally:
            csutils.xraylib.xraylib_np = keep

    @unittest.skipIf(csutils.xraylib.xraylib_np is None,
                     "xraylib_np not installed")
    def test_shape_np(self):
        self.assert_shape()

    def assert_shape(self):
        method = 'CS_Total_Kissel'
        Z, E = 1, 7
        result = csutils.eval(method, Z, E)
        self.assertTrue(isinstance(result, float))
        Z, E = 1, np.linspace(7, 7.5, 5)
        result = csutils.eval(method, Z, E)
        self.assertEqual(result.shape, (5, ))
        result, postfunc = csutils.eval(method, Z, E, applypost=False)
        self.assertEqual(result.shape, (1, 5))


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_csutils("test_shape"))
    testSuite.addTest(test_csutils("test_shape_np"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

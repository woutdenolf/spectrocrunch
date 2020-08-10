# -*- coding: utf-8 -*-

import unittest

from .. import utils
import numpy as np


class test_utils(unittest.TestCase):
    def reshape_spectrum_lines(self):
        energy = 10.0
        e, w, singlesource, singleline = utils.reshape_spectrum_lines(energy)
        self.assertTrue(singlesource)
        self.assertTrue(singleline)
        self.assertEqual(e.shape, (1, 1))
        self.assertEqual(w.shape, (1, 1))

        energy = np.array(10.0)
        e, w, singlesource, singleline = utils.reshape_spectrum_lines(energy)
        self.assertTrue(singlesource)
        self.assertTrue(singleline)
        self.assertEqual(e.shape, (1, 1))
        self.assertEqual(w.shape, (1, 1))

        energy = [10.0]
        e, w, singlesource, singleline = utils.reshape_spectrum_lines(energy)
        self.assertTrue(singlesource)
        self.assertFalse(singleline)
        self.assertEqual(e.shape, (1, 1))
        self.assertEqual(w.shape, (1, 1))

        energy = [[10.0]]
        e, w, singlesource, singleline = utils.reshape_spectrum_lines(energy)
        self.assertFalse(singlesource)
        self.assertFalse(singleline)
        self.assertEqual(e.shape, (1, 1))
        self.assertEqual(w.shape, (1, 1))


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_utils("reshape_spectrum_lines"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

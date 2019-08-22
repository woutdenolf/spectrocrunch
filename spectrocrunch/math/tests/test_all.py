# -*- coding: utf-8 -*-

import unittest
from . import test_fit2d
from . import test_fit1d
from . import test_ft
from . import test_ops
from . import test_lazy
from . import test_noisepropagation
from . import test_distributions
from . import test_interpolate


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_fit2d.test_suite())
    testSuite.addTest(test_fit1d.test_suite())
    testSuite.addTest(test_ft.test_suite())
    testSuite.addTest(test_ops.test_suite())
    testSuite.addTest(test_lazy.test_suite())
    testSuite.addTest(test_noisepropagation.test_suite())
    testSuite.addTest(test_distributions.test_suite())
    testSuite.addTest(test_interpolate.test_suite())
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

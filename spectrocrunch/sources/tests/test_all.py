# -*- coding: utf-8 -*-

import unittest
from . import test_polarization
from . import test_xray
from . import test_emspectrum


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_polarization.test_suite())
    testSuite.addTest(test_xray.test_suite())
    testSuite.addTest(test_emspectrum.test_suite())
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

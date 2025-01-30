# -*- coding: utf-8 -*-

import unittest
from . import test_base
from . import test_xray


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_base.main_test_suite())
    testSuite.addTest(test_xray.main_test_suite())
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

# -*- coding: utf-8 -*-

import unittest
from . import test_id21_ffxas
from . import test_fluoxas
from . import test_parameters


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_parameters.main_test_suite())
    testSuite.addTest(test_id21_ffxas.main_test_suite())
    testSuite.addTest(test_fluoxas.main_test_suite())
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

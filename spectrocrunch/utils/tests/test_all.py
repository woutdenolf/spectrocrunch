# -*- coding: utf-8 -*-

import unittest

from . import test_instance
from . import test_units
from . import test_classfactory
from . import test_indexing
from . import test_hashing
from . import test_signalhandling
from . import test_listtools
from . import test_lut


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_instance.main_test_suite())
    testSuite.addTest(test_units.main_test_suite())
    testSuite.addTest(test_classfactory.main_test_suite())
    testSuite.addTest(test_indexing.main_test_suite())
    testSuite.addTest(test_hashing.main_test_suite())
    testSuite.addTest(test_signalhandling.main_test_suite())
    testSuite.addTest(test_listtools.main_test_suite())
    testSuite.addTest(test_lut.main_test_suite())
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

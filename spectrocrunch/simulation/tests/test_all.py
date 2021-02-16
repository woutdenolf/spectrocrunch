# -*- coding: utf-8 -*-

import unittest

from . import test_objects
from . import test_calcnoise
from . import test_xrmc


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_objects.test_suite())
    testSuite.addTest(test_calcnoise.test_suite())
    testSuite.addTest(test_xrmc.test_suite())
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

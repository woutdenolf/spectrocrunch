# -*- coding: utf-8 -*-

import unittest
from . import test_axis
from . import test_regulargrid
from . import test_task_generic
from . import test_task_xrf

def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_axis.test_suite())
    testSuite.addTest(test_regulargrid.test_suite())
    testSuite.addTest(test_task_generic.test_suite())
    testSuite.addTest(test_task_xrf.test_suite())
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

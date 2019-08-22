# -*- coding: utf-8 -*-

import unittest
from . import test_teststack
from . import test_alignSourceDest
from . import test_transform
from . import test_align

def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_alignSourceDest.test_suite())
    testSuite.addTest(test_teststack.test_suite())
    testSuite.addTest(test_transform.test_suite())
    testSuite.addTest(test_align.test_suite())
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

import unittest
from . import test_teststack
from . import test_alignSourceDest
from . import test_transform
from . import test_align


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_alignSourceDest.main_test_suite())
    testSuite.addTest(test_teststack.main_test_suite())
    testSuite.addTest(test_transform.main_test_suite())
    testSuite.addTest(test_align.main_test_suite())
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

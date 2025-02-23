import unittest
from . import test_axis
from . import test_regulargrid
from . import test_h5merge
from . import test_task_generic
from . import test_task_xrf


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_axis.main_test_suite())
    testSuite.addTest(test_regulargrid.main_test_suite())
    testSuite.addTest(test_h5merge.main_test_suite())
    testSuite.addTest(test_task_generic.main_test_suite())
    testSuite.addTest(test_task_xrf.main_test_suite())
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

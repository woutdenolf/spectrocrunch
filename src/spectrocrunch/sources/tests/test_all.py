import unittest
from . import test_polarization
from . import test_xray
from . import test_emspectrum


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_polarization.main_test_suite())
    testSuite.addTest(test_xray.main_test_suite())
    testSuite.addTest(test_emspectrum.main_test_suite())
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

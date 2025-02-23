import unittest
from . import test_fit2d
from . import test_fit1d
from . import test_ft
from . import test_ops
from . import test_lazy
from . import test_noisepropagation
from . import test_distributions
from . import test_interpolate
from . import test_lsqlin
from . import test_quadrics
from . import test_slicing


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_fit2d.main_test_suite())
    testSuite.addTest(test_fit1d.main_test_suite())
    testSuite.addTest(test_ft.main_test_suite())
    testSuite.addTest(test_ops.main_test_suite())
    testSuite.addTest(test_lazy.main_test_suite())
    testSuite.addTest(test_noisepropagation.main_test_suite())
    testSuite.addTest(test_distributions.main_test_suite())
    testSuite.addTest(test_interpolate.main_test_suite())
    testSuite.addTest(test_lsqlin.main_test_suite())
    testSuite.addTest(test_quadrics.main_test_suite())
    testSuite.addTest(test_slicing.main_test_suite())
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

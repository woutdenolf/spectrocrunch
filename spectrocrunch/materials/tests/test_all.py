import unittest
from . import test_utils
from . import test_stoichiometry
from . import test_csutils
from . import test_compound
from . import test_mixture
from . import test_element
from . import test_visirlib
from . import test_multilayer
from . import test_pymca
from . import test_lenses
from . import test_scintillators


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_utils.main_test_suite())
    testSuite.addTest(test_stoichiometry.main_test_suite())
    testSuite.addTest(test_csutils.main_test_suite())
    testSuite.addTest(test_element.main_test_suite())
    testSuite.addTest(test_compound.main_test_suite())
    testSuite.addTest(test_mixture.main_test_suite())
    testSuite.addTest(test_visirlib.main_test_suite())
    testSuite.addTest(test_multilayer.main_test_suite())
    testSuite.addTest(test_pymca.main_test_suite())
    testSuite.addTest(test_lenses.main_test_suite())
    testSuite.addTest(test_scintillators.main_test_suite())
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

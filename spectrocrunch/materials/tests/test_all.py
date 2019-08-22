# -*- coding: utf-8 -*-

import unittest
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


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_stoichiometry.test_suite())
    testSuite.addTest(test_csutils.test_suite())
    testSuite.addTest(test_element.test_suite())
    testSuite.addTest(test_compound.test_suite())
    testSuite.addTest(test_mixture.test_suite())
    testSuite.addTest(test_visirlib.test_suite())
    testSuite.addTest(test_multilayer.test_suite())
    testSuite.addTest(test_pymca.test_suite())
    testSuite.addTest(test_lenses.test_suite())
    testSuite.addTest(test_scintillators.test_suite())
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

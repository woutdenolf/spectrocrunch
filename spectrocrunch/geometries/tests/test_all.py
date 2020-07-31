# -*- coding: utf-8 -*-

import unittest
from . import test_base
from . import test_diode
from . import test_xrf
from . import test_flatarea
from . import test_qxrf


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_base.test_suite())
    testSuite.addTest(test_diode.test_suite())
    testSuite.addTest(test_xrf.test_suite())
    testSuite.addTest(test_flatarea.test_suite())
    testSuite.addTest(test_qxrf.test_suite())
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

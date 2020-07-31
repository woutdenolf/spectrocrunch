# -*- coding: utf-8 -*-

import unittest

from . import test_localfs
from . import test_h5fs
from . import test_fs
from . import test_nxfs
from . import test_xiaedf
from . import test_spec
from . import test_excel
from . import test_h5py


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_localfs.test_suite())
    testSuite.addTest(test_h5fs.test_suite())
    testSuite.addTest(test_fs.test_suite())
    testSuite.addTest(test_nxfs.test_suite())
    testSuite.addTest(test_xiaedf.test_suite())
    testSuite.addTest(test_spec.test_suite())
    testSuite.addTest(test_excel.test_suite())
    testSuite.addTest(test_h5py.test_suite())
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

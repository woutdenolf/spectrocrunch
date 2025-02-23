import unittest

from . import test_localfs
from . import test_h5fs
from . import test_fs
from . import test_nxfs
from . import test_xiaedf
from . import test_spec
from . import test_excel
from . import test_h5py


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_localfs.main_test_suite())
    testSuite.addTest(test_h5fs.main_test_suite())
    testSuite.addTest(test_fs.main_test_suite())
    testSuite.addTest(test_nxfs.main_test_suite())
    testSuite.addTest(test_xiaedf.main_test_suite())
    testSuite.addTest(test_spec.main_test_suite())
    testSuite.addTest(test_excel.main_test_suite())
    testSuite.addTest(test_h5py.main_test_suite())
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

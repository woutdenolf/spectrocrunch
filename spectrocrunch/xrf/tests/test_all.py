# -*- coding: utf-8 -*-

import unittest


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

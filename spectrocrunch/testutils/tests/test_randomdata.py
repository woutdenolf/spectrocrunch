# -*- coding: utf-8 -*-

import unittest
import collections
import numpy as np

from .. import randomdata


class test_randomdata(unittest.TestCase):
    def test_native(self):
        for _ in range(500):
            o = randomdata.factory(types=("native",))
            # Check equality (shuffles unsorted types)
            self.assertEqual(o, o)

    def test_all(self):
        for _ in range(500):
            o = randomdata.factory()
            self.assertTrue("random" in str(o).lower())
            # Make sure raw data is generated
            data = o.data
            try:
                us = unicode(data)
            except NameError:
                us = str(data)
            except UnicodeDecodeError as e:
                us = data.decode("latin1")
            self.assertFalse("random" in us.lower())
            # Check equality (shuffles unsorted types)
            self.assertEqual(o, o)


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_randomdata("test_native"))
    testSuite.addTest(test_randomdata("test_all"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

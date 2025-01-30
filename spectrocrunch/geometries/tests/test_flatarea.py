# -*- coding: utf-8 -*-

import unittest

from .. import flatarea
from ...patch import jsonpickle


class test_flatarea(unittest.TestCase):
    def test_serialize(self):
        exclude = ()
        for name, cls in flatarea.Geometry.clsregistry.items():
            if name not in exclude:
                g1 = cls()
                g2 = jsonpickle.loads(jsonpickle.dumps(g1))
                self.assertEqual(g1, g2)


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_flatarea("test_serialize"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

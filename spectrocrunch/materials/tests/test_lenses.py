# -*- coding: utf-8 -*-

import unittest

from ..import lenses
from ...utils import units
from ...patch import jsonpickle


class test_lenses(unittest.TestCase):

    @unittest.skipIf(lenses.visirlib.PyTMM is None,
                     "PyTMM not installed")
    def test_serialize(self):
        exclude = ()
        for name, cls in lenses.Lens.clsregistry.items():
            if name not in exclude:
                l1 = cls()
                l2 = jsonpickle.loads(jsonpickle.dumps(l1))
                self.assertEqual(l1, l2)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_lenses("test_serialize"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

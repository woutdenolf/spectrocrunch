# -*- coding: utf-8 -*-

import unittest
import numpy as np

from ..import base
from ...utils import units
from ...patch import jsonpickle


class test_base(unittest.TestCase):

    def test_interpolate(self):
        o1 = base.Optics(kind='linear')
        o1.set_transmission(7, 0.2)
        o1.set_transmission(units.Quantity(7400, 'eV'), 0.8)
        def transmission(x): return o1.transmission(units.Quantity(x, 'keV'))
        self.assertEqual(transmission(7.2), 0.5)
        np.testing.assert_allclose(transmission([7.2]), [0.5])
        np.testing.assert_allclose(transmission([7.1, 7.2, 7.3]),
                                   [0.35, 0.5, 0.65])

    def test_serialize(self):
        o1 = base.Optics()
        o2 = jsonpickle.loads(jsonpickle.dumps(o1))
        self.assertEqual(o1, o2)
        o1.set_transmission([7, 7.5], [0.9, 0.85])
        o2 = jsonpickle.loads(jsonpickle.dumps(o1))
        self.assertEqual(o1, o2)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_base("test_interpolate"))
    testSuite.addTest(test_base("test_serialize"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

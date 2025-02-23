import unittest
import numpy as np

from .. import lut
from .. import units
from ...patch import jsonpickle


class test_lut(unittest.TestCase):
    def test_sort(self):
        l1 = lut.LUT()
        l1.add([3, 2], [0, -1])
        l1.add([1, 0], [-2, -3])
        self.assertEqual(l1.x.magnitude.tolist(), list(range(4)))
        self.assertEqual(l1.y.magnitude.tolist(), list(range(-3, 1)))

    def test_units(self):
        l1 = lut.LUT()
        l1.add(units.Quantity([3, 2], "cm"), units.Quantity([0, -1], "keV"))
        l1.add(units.Quantity([10, 0], "mm"), units.Quantity([-2000, -3000], "eV"))
        self.assertEqual(l1.x.magnitude.tolist(), list(range(4)))
        self.assertEqual(l1.y.magnitude.tolist(), list(range(-3, 1)))
        self.assertEqual(l1.x.units, units.ureg.Unit("cm"))
        self.assertEqual(l1.y.units, units.ureg.Unit("keV"))

    def test_interpolate(self):
        l1 = lut.LUT(kind="linear", default=units.Quantity(np.nan, "mm"))

        def func(x):
            return l1(units.Quantity(x, "keV")).to("mm").magnitude

        _ = str(l1)
        self.assertTrue(np.isnan(func(7.2)))
        self.assertTrue(np.isnan(func([7.2])).all())
        self.assertTrue(np.isnan(func([7.1, 7.2, 7.3])).all())

        l1.add(units.Quantity(7, "keV"), units.Quantity(10, "mm"))
        _ = str(l1)
        self.assertEqual(func(7.2), 10)
        np.testing.assert_allclose(func([7.2]), [10])
        np.testing.assert_allclose(func([7.1, 7.2, 7.3]), [10, 10, 10])

        l1.add(units.Quantity(7400, "eV"), units.Quantity(2, "cm"))
        _ = str(l1)
        self.assertEqual(func(7.2), 15)
        np.testing.assert_allclose(func([7.2]), [15])
        np.testing.assert_allclose(func([7.1, 7.2, 7.3]), [12.5, 15, 17.5])

    def test_serialize(self):
        l1 = lut.LUT()
        l2 = jsonpickle.loads(jsonpickle.dumps(l1))
        self.assertEqual(l1, l2)
        l1.add(1, 2)
        l2 = jsonpickle.loads(jsonpickle.dumps(l1))
        self.assertEqual(l1, l2)
        l1.add([4, 5], [6, 7])
        l2 = jsonpickle.loads(jsonpickle.dumps(l1))
        self.assertEqual(l1, l2)

    def test_add(self):
        x = units.Quantity([7, 7.1, 7.2], "keV")
        y = units.Quantity([1, 2, 3], "mA")
        l1 = lut.LUT(x=x, y=y)
        x = units.Quantity([7000, 7050, 7200], "eV")
        y = units.Quantity([4000, 5000, 6000], "uA")
        l2 = lut.LUT(x=x, y=y)
        l3 = l1 + l2

        def func(x):
            return l3(units.Quantity(x, "keV")).to("mA").magnitude.tolist()

        self.assertEqual(func([7, 7.05, 7.1, 7.2]), [5, 5, 2, 9])
        l1 += l2
        self.assertEqual(l1, l3)

    def test_replace(self):
        x = units.Quantity([7, 7.1, 7.2], "keV")
        y = units.Quantity([1, 2, 3], "mA")
        l1 = lut.LUT(x=x, y=y)
        x = units.Quantity([7000, 7050, 7200], "eV")
        y = units.Quantity([4000, 5000, 6000], "uA")
        l2 = lut.LUT(x=x, y=y)
        l1.replace(l2.x, l2.y)

        def func(x):
            return l1(units.Quantity(x, "keV")).to("mA").magnitude.tolist()

        self.assertEqual(func([7, 7.05, 7.1, 7.2]), [4, 5, 2, 6])


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_lut("test_sort"))
    testSuite.addTest(test_lut("test_units"))
    testSuite.addTest(test_lut("test_add"))
    testSuite.addTest(test_lut("test_replace"))
    testSuite.addTest(test_lut("test_interpolate"))
    testSuite.addTest(test_lut("test_serialize"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

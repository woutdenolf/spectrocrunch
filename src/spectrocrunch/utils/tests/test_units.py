import unittest
import itertools
import numpy as np

from .gendata import gendata
from .. import units
from .. import instance
from ...patch.pint import ureg


class test_units(unittest.TestCase):
    def test_quantunits(self):
        for q, arr, expand in itertools.product([True, False], repeat=3):
            # x: scalar or array
            if arr:
                x = np.arange(1.0, 10)
            else:
                x = 10.0

            # y: scalar, array, quantity(scalar), array(quantity), quantity(array)
            if q:
                unit0 = ureg.millimeter
                if arr and expand:
                    y = np.vectorize(
                        lambda a: units.Quantity(a, units=unit0), otypes=[object]
                    )(x)
                else:
                    y = units.Quantity(x, units=unit0)
                unit1 = ureg.meter
                munit = 1e-3
            else:
                unit0 = None
                unit1 = None
                munit = 1
                y = x

            # Test bubbling up and down
            if q:
                # z: quantity(scalar), quantity(array)
                z = units.Quantity(x, units=unit0)
                np.testing.assert_array_equal(units.asqarray(y), units.asqarray(z))
                # z: array(quantity)
                z = np.vectorize(
                    lambda a: units.Quantity(a, units=unit0), otypes=[object]
                )(x)
                np.testing.assert_array_equal(units.Quantity(y), units.Quantity(z))

            # Test magnitude
            a = x * munit
            b = units.magnitude(y, unit1)
            if arr:
                np.testing.assert_array_equal(a, b)
            else:
                self.assertEqual(a, b)

            # Test unit conversion
            if q:
                a = units.Quantity(x, units=unit0)
                b = units.quantity_like(y, a)
            if arr:
                np.testing.assert_array_equal(a, b)
            else:
                self.assertEqual(a, b)

    def test_asqarray(self):
        for k, v in gendata().items():
            if instance.isarray(v) or instance.isscalar(v):
                varr = units.asqarray(v)
                msg = "asqarray({}) is not an array".format(k)
                self.assertTrue(instance.isqarray(varr), msg=msg)
        for k, v in gendata().items():
            if instance.isarray(v) or instance.isscalar(v):
                barr = instance.isarray(v)
                varr, restore = units.asqarrayf(v)
                msg = "asqarrayf({}) is not an array".format(k)
                self.assertTrue(instance.isqarray(varr), msg=msg)
                v2 = restore(varr)
                barr2 = instance.isqarray(v2)
                if barr:
                    msg = "restore({}) is not an array".format(k)
                else:
                    msg = "restore({}) is an array".format(k)
                self.assertEqual(barr, barr2, msg=msg)

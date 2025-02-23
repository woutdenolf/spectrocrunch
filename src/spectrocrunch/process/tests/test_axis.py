import unittest
import numpy as np

from .. import axis
from ...utils import units


class test_axis(unittest.TestCase):
    def test_quantitative(self):
        a = units.Quantity(10, units="um")
        i = units.Quantity(13, units="um")
        j = units.Quantity(17, units="um")
        b = units.Quantity(20, units="um")
        precision = units.Quantity(1, units="nm")
        nsteps = 10

        ax1 = axis.AxisRegular(a, b, nsteps, title="x", precision=precision)
        self.assertEqual(ax1, axis.Axis(ax1.values))

        ax2 = axis.AxisRegularInc(
            a.to("mm"), ax1.stepsize, nsteps, title="y", precision=precision
        )
        self.assertEqual(ax2, axis.Axis(ax2.values))
        self.assertEqual(ax1, ax2)

        ax3 = axis.AxisNumber(a, title="z", precision=precision)
        self.assertEqual(ax3, axis.Axis(ax3.values))
        self.assertEqual(ax3, axis.AxisRegular(a, a, 0))

        ax4 = axis.factory(ax2.values, precision=ax2.precision)
        self.assertEqual(ax4, ax2)
        self.assertTrue(isinstance(ax4, axis.AxisRegular))

        o = units.Quantity(3, units="um")
        ax5 = axis.AxisSegments(
            [a, i, j, b], [nsteps, nsteps, nsteps], precision=precision
        )
        ax6 = axis.AxisSegments(
            [a - o, i, j, b], [nsteps, nsteps, nsteps], precision=precision
        )
        ax5.limits = a - o, i, j, b
        self.assertEqual(ax5, ax6)

        ax7 = axis.factory(ax5.values, precision=ax5.precision)
        self.assertEqual(ax5, ax7)
        self.assertTrue(isinstance(ax7, axis.AxisSegments))
        np.testing.assert_array_equal(ax5.limits, ax7.limits)

    def test_locate(self):
        ax = axis.factory(units.Quantity(np.arange(10) + 1, units="um"))
        self.assertEqual(ax.locate(2), 1)
        np.testing.assert_array_equal(ax.locate([2, 5]), [1, 4])
        self.assertEqual(ax.locate(-10), 0)
        ax2 = axis.factory(units.Quantity((np.arange(10) + 6) / 1000.0, units="mm"))
        self.assertEqual(ax.locate(ax2), list(range(5, 10)) + [9] * 5)

        ax = axis.factory(["list", "tuple", "set"], type="nominal")
        self.assertEqual(ax.locate("tuple"), 1)
        self.assertEqual(ax.locate(["tuple", "list"]), [1, 0])
        self.assertEqual(ax.locate("se"), 2)
        self.assertEqual(ax.locate(["se", "ist"]), [2, 0])

    def test_interpolate(self):
        ax = axis.factory(units.Quantity(np.arange(10), units="um"))
        xold, xnew = ax.interpolate(2)
        np.testing.assert_array_equal(xold, range(10))
        self.assertEqual(xnew, 2)
        xold, xnew = ax.interpolate(units.Quantity([2, 5], units="um").to("mm"))
        np.testing.assert_array_equal(xold, range(10))
        np.testing.assert_array_equal(xnew, [2, 5])

        ax = axis.factory(["a", "b", "c"], type="nominal")
        xold, xnew = ax.interpolate("b")
        self.assertEqual(xold, range(3))
        self.assertEqual(xnew, [1])
        xold, xnew = ax.interpolate(["b", "a"])
        self.assertEqual(xold, range(3))
        self.assertEqual(xnew, [1, 0])

        self.assertRaises(ValueError, ax.interpolate, "d")
        self.assertRaises(ValueError, ax.interpolate, ["a", "d"])

    def test_newlimits(self):
        x0, x1 = 0, 10
        ax1 = axis.factory(units.Quantity(np.arange(x0, x1 + 1), units="um"))
        for ai, bi in [(-2, 3), (2, -4), (3, 1), (-2, -4)]:
            a, b = x0 + ai, x1 + bi
            ax2 = axis.factory(units.Quantity(np.arange(a, b + 1), units="um"))
            self.assertEqual(ax1.newlimits(a, b), ax2)

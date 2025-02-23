import unittest

import matplotlib.pyplot as plt
import numpy as np

from .. import scene
from ...patch.pint import ureg


class test_scene(unittest.TestCase):
    def test_images(self):
        n0, n1 = 5, 10
        img = np.arange(n0 * n1).reshape(n0, n1)

        unit0 = ureg.mm
        unit1 = ureg.micrometer

        s1 = scene.Scene(unit0=unit0, unit1=unit1)

        s2 = scene.Scene(unit0=unit0, unit1=unit1)
        s2.transpose(True)
        # s2.flipx(increasing=True)
        s2.axlabels = ["dim0", "dim1"]
        s2.cmap = plt.get_cmap("gray")

        o1 = scene.Image(
            img, lim0=s1.q0([8, 8 + n0 - 1]), lim1=s1.q1([10 + n1 - 1, 10])
        )
        s1.register(o1)
        s2.register(o1)

        p0 = sorted(o1.datarange(0, border=False))
        p1 = sorted(o1.datarange(1, border=False))
        o = scene.Polyline([p0[0], p0[1], p0[1], p0[0]], [p1[0], p1[0], p1[1], p1[1]])
        s1.register(o)
        s2.register(o)
        o.set_setting("scatter", True)

        o2 = scene.Image(
            img, lim0=s1.q0([-2, -2 + n0 - 1]), lim1=s1.q1([-1, -1 + n1 - 1])
        )
        s1.register(o2)
        s2.register(o2)
        o.set_setting("scatter", True)

        p0 = sorted(o2.datarange(0, border=False))
        p1 = sorted(o2.datarange(1, border=False))
        o = scene.Text(
            [p0[0], p0[1], p0[1], p0[0]],
            [p1[0], p1[0], p1[1], p1[1]],
            labels=[1, 2, 3, 4],
        )
        s1.register(o)
        s2.register(o)

        f, ax = plt.subplots()
        s1.setaxes(ax)
        f, ax = plt.subplots()
        s2.setaxes(ax)

        # Update scene 1
        s1.updateview()

        # Shift image, axes scaling and update scene 2
        o1.lim[0] = s1.q0([9, 9 + n0 - 1])
        s2.setdatarange(0, s1.q0([0, 1]))
        s2.setdatarange(1, s1.q1([0, 1]))
        s2.updateview()
        # plt.pause(0.01)

        # Update scene 1
        s1.updateview()

        # Reset axes of scene 1
        f, ax = plt.subplots()
        s1.setaxes(ax)

        # Shift image, axes offset, different normalization and update scene 1
        o1.lim[0] = s1.q0([9, 9 + n0 - 1])

        s1.set_settings({"cnorm": "power", "cnormargs": (0.1,)})
        s1.updateview()

        # plt.pause(0.01)
        # plt.show()


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_scene("test_images"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

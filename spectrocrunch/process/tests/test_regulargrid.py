# -*- coding: utf-8 -*-

import unittest
import numpy as np

from .. import regulargrid
from .. import axis
from ...utils import units
from ...utils.tests import genindexing


class test_regulargrid(unittest.TestCase):
    def _generate(self, stackaxis, shape=(2, 3, 4), stackdim=0):
        dtype = np.float32
        stackaxis = axis.factory(stackaxis, type="ordinal")
        axes = [axis.factory(range(n)) for n in shape]
        axes.insert(stackdim, stackaxis)
        shape = [len(ax) for ax in axes]
        data = np.random.normal(size=shape).astype(dtype) * 100
        return regulargrid.RegularGrid(axes, data, stackdim=stackdim)

    def test_indexing(self):
        stackaxis = ["Fe-K", "Si-K", "Al-K", "S-K", "Ce-L"]
        for i in range(3):
            grid = self._generate(stackaxis, shape=(6, 8), stackdim=i)
            self._check_grid(grid)

    def test_interpolate(self):
        for degree in [0, 1, 2]:
            for asgrid in [True, False]:
                rtol = 1e-3

                # 1D
                stackaxis = ["Fe-K", "Si-K", "Al-K", "S-K", "Ce-L"]
                grid = self._generate(stackaxis, shape=(), stackdim=0)
                data = grid.interpolate("Fe-K", degree=degree, asgrid=asgrid)
                np.testing.assert_allclose(grid.data[0], data, rtol=rtol)
                data = grid.interpolate(["Fe-K", "Ce-L"], degree=degree, asgrid=asgrid)
                np.testing.assert_allclose(grid.data[[0, -1]], data, rtol=rtol)

                # 2D
                grid = self._generate(stackaxis, shape=(6,), stackdim=0)
                data = grid.interpolate("Fe-K", None, degree=degree)
                np.testing.assert_allclose(grid.data[np.newaxis, 0, :], data, rtol=rtol)
                data = grid.interpolate(
                    ["Fe-K", "S-K", "Ce-L"], [0, 2, 3], degree=degree, asgrid=asgrid
                )
                ind = [0, 3, 4], [0, 2, 3]
                if asgrid:
                    ind = tuple(np.meshgrid(*ind, indexing="ij"))
                np.testing.assert_allclose(grid.data[ind], data, rtol=rtol)

                # 3D
                grid = self._generate(stackaxis, shape=(6, 8), stackdim=1)
                data = grid.interpolate(None, "Fe-K", None, degree=degree)
                np.testing.assert_allclose(
                    grid.data[:, np.newaxis, 0, :], data, rtol=rtol
                )
                data = grid.interpolate(
                    [1, 3, 4],
                    ["Fe-K", "S-K", "Ce-L"],
                    [0, 2, 3],
                    degree=degree,
                    asgrid=asgrid,
                )
                ind = [1, 3, 4], [0, 3, 4], [0, 2, 3]
                if asgrid:
                    ind = tuple(np.meshgrid(*ind, indexing="ij"))
                np.testing.assert_allclose(grid.data[ind], data, rtol=rtol)

    def _check_grid(self, grid):
        data = grid.values
        self.assertEqual(grid.shape, data.shape)
        self.assertEqual(grid.ndim, data.ndim)
        self.assertEqual(grid.size, data.size)
        np.testing.assert_array_equal(grid[:], data)

        indices = genindexing.genindexingn(
            data.shape, advanced=False, eco=False, nmax=50
        )
        for index in indices:
            np.testing.assert_array_equal(grid[index], data[index])

        for a, b in zip(grid, data):
            np.testing.assert_array_equal(a, b)


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_regulargrid("test_indexing"))
    testSuite.addTest(test_regulargrid("test_interpolate"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

# -*- coding: utf-8 -*-

import unittest
from spectrocrunch.process import h5merge
import h5py
import numpy
import numpy.testing


class test_h5merge(unittest.TestCase):
    def test_tile(self):
        scan_shape = (7, 9)
        tile_shape = (2, 2)
        shapes = [(63,), (62,), (61,)]
        shape_map = {s: scan_shape for s in shapes}

        expected = numpy.zeros((14, 18), dtype=int)

        off = 1
        expected[0:7, 0:9] = numpy.arange(63).reshape((7, 9)) + off
        off += 63
        expected[0:6, 9:18] = numpy.arange(54).reshape((6, 9)) + off
        off += 62
        expected[7:13, 0:9] = numpy.arange(54).reshape((6, 9)) + off

        with h5py.File("test.h5", "w") as f:
            sources = list()
            off = 1
            for i, shape in enumerate(shapes):
                group = f.create_group(f"group{i}")
                sources.append(group)
                n = numpy.prod(shape, dtype=int)
                data = numpy.arange(off, off + n)
                off += n
                group.create_dataset("data", data=data.reshape(shape))

            merged = f.create_group("merged")

            h5merge.merge_h5groups(
                merged,
                "group0-2",
                sources,
                shape_map,
                nscandim=len(scan_shape),
                tile_shape=tile_shape,
            )
            vds = merged["group0-2"]["data"][()]
            numpy.testing.assert_array_almost_equal(expected, vds)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_h5merge("test_tile"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

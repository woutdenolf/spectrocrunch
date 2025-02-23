import os
import unittest

import h5py
import numpy
import numpy.testing
from testfixtures import TempDirectory

from spectrocrunch.process import h5merge


class test_h5merge(unittest.TestCase):
    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()

    def test_tile(self):
        scan_shape = (9, 7)  # F-order
        tile_shape = (2, 3)  # F-order

        # only 5 of 6 maps, 2 incomplete
        shapes = [(63,), (62,), (61,), (63,), (63,)]
        shape_map = {s: scan_shape for s in shapes}

        expected = numpy.zeros((21, 18), dtype=int)

        off = 1
        expected[0:7, 0:9] = numpy.arange(63).reshape((7, 9)) + off
        off += 63
        expected[0:6, 9:18] = numpy.arange(54).reshape((6, 9)) + off
        off += 62
        expected[7:13, 0:9] = numpy.arange(54).reshape((6, 9)) + off
        off += 61
        expected[7:14, 9:18] = numpy.arange(63).reshape((7, 9)) + off
        off += 63
        expected[14:21, 0:9] = numpy.arange(63).reshape((7, 9)) + off

        with h5py.File(os.path.join(self.dir.path, "test.h5"), "w") as f:
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
            print(expected)
            print()
            print(vds)
            numpy.testing.assert_array_almost_equal(expected, vds)

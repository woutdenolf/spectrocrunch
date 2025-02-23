import unittest

import numpy
from .. import slicing


class test_slicing(unittest.TestCase):
    def test_slice_generator(self):
        itemsize = 4  # bytes
        mb_threshold = 100  # MB
        n_threshold = mb_threshold * 1024**2 // itemsize

        shape = (n_threshold + 1000, 2, 2, 4)
        nslices = len(
            list(
                slicing.slice_generator(shape, numpy.float32, mb_threshold=mb_threshold)
            )
        )
        self.assertEqual(nslices, 2 * 2 * 4)

        p = n_threshold**0.5
        shape = (p, p - 1000, 2, 4)
        nslices = len(
            list(
                slicing.slice_generator(shape, numpy.float32, mb_threshold=mb_threshold)
            )
        )
        self.assertEqual(nslices, 2 * 4)

        p = n_threshold ** (1 / 3.0)
        shape = (p, p, p - 100, 4)
        nslices = len(
            list(
                slicing.slice_generator(shape, numpy.float32, mb_threshold=mb_threshold)
            )
        )
        self.assertEqual(nslices, 4)

        p = n_threshold ** (1 / 4.0)
        shape = (p, p, p, p - 10)
        nslices = len(
            list(
                slicing.slice_generator(shape, numpy.float32, mb_threshold=mb_threshold)
            )
        )
        self.assertEqual(nslices, 1)

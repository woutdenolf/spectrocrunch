# -*- coding: utf-8 -*-

import itertools
import numpy


def calc_nbytes(dtype):
    return numpy.array(0, dtype=dtype).nbytes


def find_slice_axis(shape, dtype, mb_threshold=100):
    mb = numpy.cumprod(shape) * calc_nbytes(dtype) / 1024 ** 2
    axis = numpy.where(mb <= mb_threshold)[0]
    if axis.size:
        return axis[-1]
    else:
        return 0


def slice_generator(shape, dtype, mb_threshold=100):
    axis = find_slice_axis(shape, dtype, mb_threshold=mb_threshold)
    choices = [[slice(None)] for _ in range(axis + 1)] + [
        range(shape[i]) for i in range(axis + 1, len(shape))
    ]
    for idx in itertools.product(*choices):
        yield idx

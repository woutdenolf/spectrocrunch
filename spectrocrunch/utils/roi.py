# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import operator


def cliproi(shape, roi):
    """Make sure that a ROI does not exceeds the maximal size.

    Args:
        shape (n-tuple): array shape (n1, n2, ...)
        roi (n-2-tuple): array range indices ((a1,b1),(a2,b2),...)

    Returns:
        n-2-list: clipped ROI [[a1,b1],[a2,b2],...]
    """
    if len(shape) != len(roi):
        raise ValueError("Dimensions for shape and ROI should be the same")

    roinew = []

    for n, (a, b) in zip(shape, roi):
        if a is None:
            a = 0
        else:
            if a < 0:
                a += n
            a = max(0, min(a, n-1))

        if b is None:
            b = n
        else:
            if b < 0:
                b += n
            b = max(0, min(b, n))

        roinew += [[a, b]]

    return roinew


def mergeroi1d(intervals):
    """
    Args:
        intervals(list(2-tuple)):

    Returns:
        generator
    """
    sorted_intervals = sorted(intervals, key=operator.itemgetter(0))
    if not sorted_intervals:
        return

    low, high = sorted_intervals[0]
    for iv in sorted_intervals[1:]:
        if iv[0] <= high:
            high = max(high, iv[1])
        else:
            yield low, high
            low, high = iv

    yield low, high

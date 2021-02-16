# -*- coding: utf-8 -*-

import numpy


def cramer(A, b):
    # A.x = b
    detA = numpy.linalg.det(A)
    if detA == 0:
        raise RuntimeError("Singular matrix")
    nrow, ncol = A.shape

    sol = [None] * ncol
    for i in range(ncol):
        Ai = numpy.array(A, copy=True)
        Ai[:, i] = b
        sol[i] = numpy.linalg.det(Ai)

    return numpy.asarray(sol) / detA

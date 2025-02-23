import random
import itertools
import numpy as np
import numbers

from .. import indexing


def genindexing(dim, advanced=False, eco=False):

    if dim > 1:
        a = random.randint(-dim, dim - 1)
        i = random.randint(1, dim - 1)
        j = random.randint(-dim, -1)
        k = random.randint(1, max(1, dim // 2))
    else:
        a = 0
        i = 0
        j = 0
        k = 1

    if eco:
        ret = [a, slice(i, j, k)]
    else:
        ret = [
            0,
            i,
            j,
            slice(i, j, k),
            slice(None, j, k),
            slice(i, None, k),
            slice(i, j, None),
            slice(None, None, k),
            slice(i, None, None),
            slice(None, j, None),
            slice(None),
            slice(0),
        ]

    # Advanced indexing
    if advanced:
        n = min(dim // 2, 2)
        iarr = [random.randint(-dim, dim - 1) for i in range(n)]
        barr = [True] * n + [False] * (n + n % 2)
        random.shuffle(iarr)
        random.shuffle(barr)

        if eco:
            ret += [iarr, barr, Ellipsis, np.newaxis]
        else:
            ret += [[], [a], iarr, barr, Ellipsis, np.newaxis]
    return ret


def valid(index, shape):
    # No more than one '...'
    if sum([i is Ellipsis for i in index]) > 1:
        return False

    # More than one list: they must have the same length
    b = [i for i in index if isinstance(i, list)]
    if len(b) > 1:
        b = [len(i) for i in b]
        if b.count(b[0]) != len(b):
            return False

    # Check dimensions
    index = indexing.expand(index, len(shape))
    indexnonone = [ind for ind in index if ind is not None]
    b = True
    for i, ind in enumerate(indexnonone):
        if isinstance(ind, list):
            if len(ind) == 0:
                continue

            # Boolean array must have the correct size
            if all(isinstance(j, bool) for j in ind) and len(ind) > 0:
                b &= len(ind) == shape[i]

            # Integer array musn't go beyond boundaries
            if all(isinstance(j, numbers.Number) for j in ind):
                b &= max(ind) < shape[i] and min(ind) >= -shape[i]

        elif isinstance(ind, numbers.Number):
            b &= ind < shape[i] and ind >= -shape[i]

    return b


def genindexingn(shape, advanced=False, nmax=None, eco=False):
    lst = [genindexing(s, advanced=advanced, eco=eco) for s in shape]
    lst = list(itertools.product(*lst))
    if nmax is not None:
        if len(lst) > nmax:
            random.shuffle(lst)
            lst = lst[:nmax]
    return [t for t in lst if valid(t, shape)]

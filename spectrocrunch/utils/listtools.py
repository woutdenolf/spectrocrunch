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

import collections
import operator
import itertools
import numpy as np

from . import instance


def flatten(l):
    """Flatten iterables

    Args:
        l(anything):
    Returns:
        list
    """
    if instance.isiterable(l) and not instance.isstring(l):
        try:
            it = iter(l)
        except TypeError:
            yield l
        else:
            for el in it:
                for el2 in flatten(el):
                    yield el2
    else:
        yield l


def numpy_flatten(l):
    return np.asarray(list(flatten(l)))


def listadvanced_bool(lst, barr, bnot=False):
    """Advanced list indexing: boolean array

    Args:
        lst(list):
        barr(array or bool): array of booleans
    Returns:
        list
    """
    if bnot:
        barr = map(operator.not_, barr)
    return list(itertools.compress(lst, barr))


def listadvanced_int(lst, ind):
    """Advanced list indexing: integer array

    Args:
        lst(list):
        ind(array):
    Returns:
        list
    """
    return [lst[i] for i in ind]


def listadvanced(lst, ind):
    """Advanced list indexing: integer or bool array

    Args:
        lst(list):
        ind(array):
    Returns:
        list
    """
    if instance.isboolsequence(ind):
        return listadvanced_bool(lst, ind)
    else:
        return listadvanced_int(lst, ind)


def where(lst, func):
    """Indices are particular elements

    Args:
        lst(list):
        func(callable): one argument
    Returns:
        list
    """
    return [i for i, l in enumerate(lst) if func(l)]


def sort2lists(list1, list2):
    """Sort list1 and list2 based on list1

    Args:
        list1(list):
        list2(list):
    Returns:
        list,list
    """
    return tuple(list(t) for t in zip(*sorted(zip(list1, list2),
                                              key=operator.itemgetter(0))))


def unique2lists(list1, list2, add=False):
    """Unique list1 and list2 based on list1
    Args:
        list1(list):
        list2(list):
        add(Optional(bool)): add list2 elements with list1 duplicates
    Returns:
        list,list
    """
    if add:
        cntr = collections.Counter()

        def cntr_add(x, y):
            b = x not in cntr
            cntr[x] += y
            return b
        list1 = [x1 for x1, x2 in zip(list1, list2) if cntr_add(x1, x2)]
        list2 = [cntr[x1] for x1 in list1]
        return list1, list2
    else:
        seen = set()
        seen_add = seen.add
        list1, list2 = tuple(zip(*[[x1, x2] for x1, x2 in zip(list1, list2)
                                   if not (x1 in seen or seen_add(x1))]))
        return list(list1), list(list2)


def sumrepeats(labels, counts):
    """

    Args:
        list1(list):
        list2(list):
    Returns:
        list,list
    """
    c = collections.Counter()
    for l, cnt in zip(labels, counts):
        c.update({l: cnt})
    return c.keys(), c.values()


def swap(lst, i, j):
    if i != j:
        lst[i], lst[j] = lst[j], lst[i]
    return lst


def roll(lst, n):
    if n != 0:
        n = abs(n)
        lst = list(itertools.islice(itertools.cycle(lst), n, n+len(lst)))
    return lst


def move(lst, i, j):
    if i != j:
        lst.insert(j, lst.pop(i))
    return lst


def length(x):
    try:
        return len(x)
    except TypeError:
        return 1


def aslist(x):
    try:
        return list(x)
    except:
        return [x]


def filterfalse(predicate, iterable):
    # filterfalse(lambda x: x%2, range(10)) --> 0 2 4 6 8
    if predicate is None:
        predicate = bool
    for x in iterable:
        if not predicate(x):
            yield x


def unique_everseen(iterable, key=None):
    "List unique elements, preserving order. Remember all elements ever seen."
    # unique_everseen('AAAABBBCCDAABBB') --> A B C D
    # unique_everseen('ABBCcAD', str.lower) --> A B C D
    seen = set()
    seen_add = seen.add
    if key is None:
        for element in filterfalse(seen.__contains__, iterable):
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element

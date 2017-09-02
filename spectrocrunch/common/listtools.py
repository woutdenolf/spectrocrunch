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

from . import instance



def flatten(l):
    """Flatten list

    Args:
        l(list):
    Returns:
        list
    """
    for el in l:
        if isinstance(el, collections.Iterable) and not instance.isstring(el):
            for sub in flatten(el):
                yield sub
        else:
            yield el

def listadvanced_bool(lst,barr,bnot=False):
    """Advanced list indexing: boolean array

    Args:
        lst(list):
        barr(array or bool): array of booleans
    Returns:
        list
    """
    if len(lst)!=len(barr):
        raise IndexError("boolean index did not match indexed list; length is {} but boolean dimension is {}".format(len(lst),len(barr)))
    if bnot:
        return [item for b,item in zip(barr,lst) if not b]
    else:
        return [item for b,item in zip(barr,lst) if b]

def listadvanced_int(lst,ind):
    """Advanced list indexing: integer array

    Args:
        lst(list):
        ind(array):
    Returns:
        list
    """
    return [lst[i] for i in ind]

def listadvanced(lst,ind):
    """Advanced list indexing: integer or bool array

    Args:
        lst(list):
        ind(array):
    Returns:
        list
    """
    if instance.isboollist(ind):
        return listadvanced_bool(lst,ind)
    else:
        return listadvanced_int(lst,ind)

def where(lst,func):
    """Indices are particular elements

    Args:
        lst(list):
        func(callable): one argument
    Returns:
        list
    """
    return [i for i,l in enumerate(lst) if func(l)]
    
def sort2lists(list1, list2):
    """Sort list1 and list2 based on list1

    Args:
        list1(list):
        list2(list):
    Returns:
        list,list
    """
    return tuple(list(t) for t in itertools.izip( *sorted(itertools.izip(list1, list2),key=operator.itemgetter(0)) ))

def weightedsum(labels, counts):
    """

    Args:
        list1(list):
        list2(list):
    Returns:
        list,list
    """
    c = collections.Counter()
    for l,cnt in itertools.izip(labels,counts):
        c.update({l:cnt})
    return c.keys(),c.values()
        

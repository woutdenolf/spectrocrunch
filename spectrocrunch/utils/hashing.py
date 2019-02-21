# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
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

import hashlib
import json
import itertools
from ..patch import jsonpickle
from ..utils import instance
from ..utils import listtools


def calchash(x):
    """
    Deterministic MD5 hash of data:
        - handle unsorted data types recursively (dict, set, frozenset)
        - 

    Returns:
        str
    """
    #print('calchash({})'.format(x))
    if x:
        # Convert x to a list of hashes
        if instance.isarray(x):
            if instance.isset(x):
                x = sorted(calchash(y) for y in x)
            else:
                x = [calchash(y) for y in x]
        elif instance.ismapping(x):
            #print('keys before: {}'.format(x.keys()))
            keys = [calchash(k) for k in x.keys()]
            values = [calchash(v) for v in x.values()]
            if not instance.isorderedmapping(x):
                keys, values = listtools.sort2lists(keys, values)
            #print('keys after: {}'.format(keys))
            x = keys + values
        elif hasattr(x, '__getstate__'):
            x = [calchash(x.__getstate__())]
        # Convert x to a unique string
        #print('jsonpickle: {}'.format(x))
        x = jsonpickle.encode(x)  # str
    else:
        #print('type: {}'.format(x))
        x = 'type().__name__==' + type(x).__name__
    # Convert unicode to bytes if needed
    if not isinstance(x, bytes):
        x = x.encode()
    # MD5 hash (str) of x(bytes)
    return hashlib.md5(x).hexdigest()


def hashequal(a, b):
    return calchash(a) == calchash(b)


def mergehash(*hashes):
    if len(hashes) > 1:
        return calchash(hashes)
    else:
        return hashes[0]

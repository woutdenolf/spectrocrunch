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
import itertools
from ..patch import jsonpickle
from ..utils import instance
from ..utils import listtools


def string_repr_to_hash(x):
    if not isinstance(x, bytes):
        return x.encode('utf-8')
    return x


def classname_repr_to_hash(classname):
    return string_repr_to_hash("<class '{}'>".format(classname))


def class_repr_to_hash(o):
    return classname_repr_to_hash(o.__class__.__name__)


def any_repr_to_hash(x):
    return string_repr_to_hash(jsonpickle.dumps(x, sort_keys=True))


def getstate(x):
    return jsonpickle.flatten(x)


def calchash(x, _depth=0):
    """
    Deterministic MD5 hash of data:

    - handle unsorted data types recursively (dict, set, frozenset)
    - ignore string type difference when possible
    - ignore number type difference when possible (True == 1 == 1. == 1+0j)
    - uses jsonpickle to get the state of custom objects

    Returns:
        bytes: can be 'ascii' decoded
    """
    if _depth == 100:
        raise RuntimeError('Maximal hash recursion depth is reached')
    else:
        _depth += 1
    # Convert x to a list(bytes)
    if x is None:
        x = [class_repr_to_hash(x)]
    elif instance.isstring(x):
        x = [classname_repr_to_hash('string'), string_repr_to_hash(x)]
    elif instance.isnumber(x):
        try:
            intx = int(x)
            if intx == x:
                x = intx
        except:
            pass
        x = [classname_repr_to_hash('number'), string_repr_to_hash(str(x))]
    elif instance.isarray(x):
        if instance.isarray0(x):
            x = [calchash(x.item(), _depth=_depth)]
        elif len(x):
            if instance.isset(x):
                x = sorted(calchash(y, _depth=_depth) for y in x)
            else:
                x = [calchash(y, _depth=_depth) for y in x]
        else:
            # Empty array
            x = [class_repr_to_hash(x)]
    elif instance.ismapping(x):
        if len(x):
            keys = [calchash(k, _depth=_depth) for k in x.keys()]
            values = [calchash(v, _depth=_depth) for v in x.values()]
            if not instance.isorderedmapping(x):
                keys, values = listtools.sort2lists(keys, values)
            x = keys + values
        else:
            # Empty mapping
            x = [class_repr_to_hash(x)]
    elif instance.isquantity(x):
        x = [calchash(x.magnitude, _depth=_depth),
             calchash(str(x.units), _depth=_depth)]
    else:
        x = [calchash(getstate(x), _depth=_depth)]
    # MD5 hash of list of bytes
    #print(x)
    return string_repr_to_hash(hashlib.md5(b''.join(x)).hexdigest())


def hashequal(a, b):
    return calchash(a) == calchash(b)


def mergehash(*hashes):
    if len(hashes) > 1:
        return calchash(hashes)
    else:
        return hashes[0]

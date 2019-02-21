# -*- coding: utf-8 -*-
import collections
from array import array
import numpy as np
from ...patch.pint import ureg


def gendata():
    data = {}
    data['type'] = type
    data['none'] = None
    data['bool'] = True
    data['int'] = int(1)
    data['float'] = float(1)
    data['str'] = 'abc'
    data['unicode'] = u'äbc'
    tmp = b'abc'
    data['bytes'] = tmp
    data['bytearray'] = bytearray(tmp)
    data['array'] = array('b', tmp)
    tmp = 1, 2
    data['tuple'] = tuple(tmp)
    data['list'] = list(tmp)
    try:
        data['xrange'] = xrange(5)
    except NameError:
        data['xrange'] = range(5)
    data['frozenset'] = frozenset(tmp)
    data['set'] = set(tmp)
    data['ndarray'] = np.zeros(5)
    data['ndarray0'] = np.array(5)
    data['npscalar'] = np.float32(5)
    class Dummy:
        pass
    data['nparray'] = np.array([Dummy()]*5)
    data['nparray0'] = np.array(Dummy())
    data['qarray'] = ureg.Quantity(np.zeros(5), ureg.meter)
    data['qarray0'] = ureg.Quantity(np.array(5), ureg.meter)
    data['qscalar'] = ureg.Quantity(np.float32(5), ureg.meter)
    tmp = [('yellow', 1), ('blue', 2)]
    data['dict'] = dict(tmp)
    data['odict'] = collections.OrderedDict(tmp)
    data['ddict'] = collections.defaultdict(tuple, tmp)
    return data

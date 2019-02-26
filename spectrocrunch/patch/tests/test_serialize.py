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

import unittest
import os
import numpy as np
import uncertainties
from uncertainties import umath

from ...utils import units
from .. import jsonpickle


def equal(a, b):
    if type(a) == np.ndarray:
        return np.array_equal(a, b)
    else:
        return a == b


class ExampleClass(object):

    def __init__(self, x, y, z=None):
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self, other):
        return equal(self.x, other.x) and \
            equal(self.y, other.y) and \
            equal(self.z, other.z)

    def __getstate__(self):
        return self.__dict__.copy()

    def __setstate__(self, state):
        self.__dict__.update(state)


class LinkClass1(object):

    def __init__(self, *links):
        self.links = links

    def __repr__(self):
        return '{}->{}'.format(
            id(self),
            repr(self.links)
        )


class LinkClass2(LinkClass1):
    pass


class LinkClassHandler(jsonpickle.BaseHandler):

    def flatten(self, obj, state):
        state['links'] = self.context.flatten(obj.links, reset=False)
        return state

    def restore(self, state):
        links = state['links']
        links = self.context.restore(links, reset=False)
        return LinkClass2(links)


jsonpickle.register(LinkClass2, LinkClassHandler)


class test_serialize(unittest.TestCase):

    def test_units(self):
        a = units.Quantity(5, "keV")
        b = 10.
        c = ExampleClass(10, [10, 20], z=30)
        d = ExampleClass(np.array([1, 2, 3]),
                         ExampleClass(-10, -20, z=-30),
                         z=ExampleClass(-10, ExampleClass(-10, -20, z=-30),
                                        z=ExampleClass(-10, -20, z=-30)))
        data = {"a": a, "b": b, "c": c, "d": d}
        serialized = jsonpickle.encode(data)
        deserialized = jsonpickle.decode(serialized)
        self.assertEqual(data, deserialized)

    @unittest.skip('cyclic jsonpickle bug')
    def test_uncertainties(self):
        data = uncertainties.ufloat(1, 2)
        serialized = jsonpickle.encode(data)
        deserialized = jsonpickle.decode(serialized)
        print(serialized)
        print(deserialized)
        self.assertEqual(data.nominal_value, deserialized.nominal_value)
        self.assertEqual(data.std_dev, deserialized.std_dev)
        data += 2
        serialized = jsonpickle.encode(data)
        deserialized = jsonpickle.decode(serialized)
        print(serialized)
        print(deserialized)
        self.assertEqual(data.nominal_value, deserialized.nominal_value)
        self.assertEqual(data.std_dev, deserialized.std_dev)

    @unittest.skip('cyclic jsonpickle bug')
    def test_cyclic(self):
        for _class in (LinkClass2, ):
            child = _class(None)
            instance = _class(child, child)
            child.links = instance,
            print(instance)
            self.assertEqual(id(instance.links[0]), id(instance.links[1]))
            self.assertEqual(id(instance.links[0].links[0]), id(instance))
            print(jsonpickle.Pickler().flatten(instance))
            encoded = jsonpickle.encode(instance)
            instance = jsonpickle.decode(encoded)
            print(instance)
            self.assertEqual(id(instance.links[0]), id(instance.links[1]))
            self.assertEqual(id(instance.links[0].links[0]), id(instance))
            continue

            a = {}
            b = LinkClass1(a)
            c = LinkClass1(b)
            a['link'] = b
            self.assertEqual(id(c.links[0]), id(b))
            self.assertEqual(id(b.links[0]), id(a))
            print(jsonpickle.Pickler().flatten(c))
            serialized = jsonpickle.encode(c)
            c = jsonpickle.decode(serialized)
            b = c.links[0]
            a = b.links[0]
            self.assertEqual(id(c.links[0]), id(b))
            self.assertEqual(id(b.links[0]), id(a))


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_serialize("test_cyclic"))
    testSuite.addTest(test_serialize("test_units"))
    testSuite.addTest(test_serialize("test_uncertainties"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

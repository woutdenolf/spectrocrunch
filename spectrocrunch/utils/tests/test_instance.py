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
import numpy as np

from .gendata import gendata
from .. import units
from .. import instance
from ...math import noisepropagation


class test_instance(unittest.TestCase):

    def test_instance(self):
        select = 'str', 'unicode', 'bytes', 'npunicode', 'npbytes', 'npstr'
        self._check_instance(instance.isstring, select)
        select = 'npunicode', 'npbytes', 'npstr'
        self._check_instance(instance.isnpstring, select)
        select = ('str_list', 'unicode_list', 'bytes_list',
                  'npstr_array', 'npunicode_array', 'npbytes_array')
        self._check_instance(instance.isstringarray, select)
        select = ('tuple', 'list', 'xrange',
                  'str_list', 'unicode_list', 'bytes_list')
        self._check_instance(instance.issequence, select)
        select = 'set', 'frozenset'
        self._check_instance(instance.isset, select)
        select = ('tuple', 'list', 'xrange',
                  'set', 'frozenset',
                  'array', 'bytearray',
                  'ndarray', 'ndarray0',
                  'nparray', 'nparray0',
                  'qarray', 'qarray0',
                  'uarray', 'uarray0',
                  'str_list', 'unicode_list', 'bytes_list',
                  'npstr_array', 'npunicode_array', 'npbytes_array',
                  'dict', 'odict', 'ddict')
        self._check_instance(instance.isiterable, select)
        select = ('tuple', 'list', 'xrange',
                  'set', 'frozenset',
                  'array', 'bytearray',
                  'ndarray', 'ndarray0',
                  'nparray', 'nparray0',
                  'qarray', 'qarray0',
                  'uarray', 'uarray0',
                  'str_list', 'unicode_list', 'bytes_list',
                  'npstr_array', 'npunicode_array', 'npbytes_array')
        self._check_instance(instance.isarray, select)
        select = 'ndarray0', 'nparray0', 'qarray0', 'uarray0'
        self._check_instance(instance.isarray0, select)
        select = ('tuple', 'list', 'xrange',
                  'set', 'frozenset',
                  'array', 'bytearray',
                  'ndarray', 'nparray',
                  'qarray', 'uarray',
                  'str_list', 'unicode_list', 'bytes_list',
                  'npstr_array', 'npunicode_array', 'npbytes_array')
        self._check_instance(instance.isarraynot0, select)
        select = ('bool', 'int', 'float',
                  'nan', 'inf', 'npint', 'npfloat',
                  'npunicode', 'npbytes', 'npstr',
                  'qscalar', 'uscalar')
        self._check_instance(instance.isscalar, select)
        select = 'qscalar',
        self._check_instance(instance.isqscalar, select)
        select = ('npint', 'npfloat',
                  'npunicode', 'npbytes', 'npstr')
        self._check_instance(instance.isnpscalar, select)
        select = 'uscalar',
        self._check_instance(instance.isuscalar, select)
        select = ('bool', 'int', 'float',
                  'nan', 'inf', 'npint', 'npfloat')
        self._check_instance(instance.isnumber, select)
        select = ('int', 'float',
                  'nan', 'inf', 'npint', 'npfloat')
        self._check_instance(instance.isnumeric, select)
        select = 'npint', 'npfloat'
        self._check_instance(instance.isnpnumber, select)
        select = 'bool', 'int', 'npint'
        self._check_instance(instance.isinteger, select)
        select = 'qarray', 'qarray0', 'qscalar'
        self._check_instance(instance.isquantity, select)
        select = ('ndarray', 'ndarray0',
                  'nparray', 'nparray0',
                  'uarray', 'uarray0',
                  'npstr_array', 'npunicode_array', 'npbytes_array')
        self._check_instance(instance.isnparray, select)
        select = ('ndarray', 'ndarray0',
                  'nparray', 'nparray0',
                  'uarray', 'uarray0',
                  'npint', 'npfloat',
                  'npstr', 'npunicode', 'npbytes',
                  'npstr_array', 'npunicode_array', 'npbytes_array')
        self._check_instance(instance.isnumpy, select)
        select = 'qarray', 'qarray0'
        self._check_instance(instance.isqarray, select)
        select = 'uarray', 'uarray0'
        self._check_instance(instance.isuarray, select)
        select = 'dict', 'ddict', 'odict'
        self._check_instance(instance.ismapping, select)
        select = 'odict',
        self._check_instance(instance.isorderedmapping, select)

    def _check_instance(self, isinstance, select):
        for k, v in gendata().items():
            b1 = isinstance(v)
            b2 = k in select
            msg = '{}({}) != {}'.format(isinstance.__name__, k, b2)
            self.assertEqual(b1, b2, msg=msg)

    def _test_asarray(self, x):
        if instance.isarray(x) and not instance.isquantity(x):
            xtype = np.ndarray
        else:
            xtype = type(x)
        y, func = instance.asarrayf(x)
        self.assertTrue(isinstance(y, np.ndarray))
        y = func(y)
        self.assertEqual(type(y), xtype)

    def test_asarray(self):
        for k, v in gendata().items():
            varr = instance.asarray(v)
            self.assertTrue(instance.isarray(varr))
        for k, v in gendata().items():
            barr = instance.isarray(v)
            varr, restore = instance.asarrayf(v)
            barr2 = instance.isarray(restore(varr))
            self.assertEqual(barr, barr2)

        a = np.int64(0)
        b = np.int64(1)
        nums = [a, np.array(a), [a], np.array([a]), [a, b], np.array([a, b])]
        checktypes = [True, True, False, True, False]
        for num, check in zip(nums, checktypes):
            check = True
            x = num
            self._test_asarray(x)
            x = units.Quantity(num, units=units.ureg.millimeter)
            self._test_asarray(x)
            x = np.vectorize(lambda a: units.Quantity(
                a, units=units.ureg.millimeter), otypes=[object])(x)
            self._test_asarray(x)
            x = noisepropagation.randomvariable(num, num)
            self._test_asarray(x)

        objs = ["abc", np.array("abc"), np.array(
            ["abc"]), ["abc"], ("abc", 1, 2)]
        checktypes = [True, True, True, False, False]
        for o, check in zip(nums, checktypes):
            self._test_asarray(x)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_instance("test_instance"))
    testSuite.addTest(test_instance("test_asarray"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

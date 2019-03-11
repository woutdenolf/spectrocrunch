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
from .. import hashing
from ...patch.pint import ureg
from ...testutils.randomdata import factory
from .. import instance


class test_hashing(unittest.TestCase):

    def test_random(self):
        for _ in range(500):
            #print('+++++++++++++++')
            o = factory()
            a = o.data
            b = o.data
            try:
                ha = hashing.calchash(a)
            except UnicodeDecodeError:
                print(a)
                raise
            try:
                hb = hashing.calchash(b)
            except UnicodeDecodeError:
                print(b)
                raise
            if ha != hb:
                print('---')
                print(a)
                print(b)
                assert()
            self.assertTrue(hashing.hashequal(a, b))

    def test_string(self):
        # Empty string
        hash = b'7d52e4059c583de024b90a033ef4577e'
        self._assert_hash(b'', hash)
        self._assert_hash(u'', hash)
        # Ascii
        hash = b'4a77c5fa0df8eeb6f66ea3e95148e57e'
        self._assert_hash(b'abc', hash)
        self._assert_hash(u'abc', hash)
        # Extended ascii
        hash = b'9af537523552e3988d1eadc7dd4eb048'
        self._assert_hash(b'\xe423', hash)
        # UTF8
        hash = b'3a00931200bada8bb1a4c4c36185898f'
        self._assert_hash(u'\u0101bc', hash)

    def test_number(self):
        hash = b'cbeddc349fedb70a755c2b345e0f6a93'
        self._assert_hash(float('nan'), hash)
        hash = b'a713229c834747840d40cafa239c52ba'
        self._assert_hash(float('inf'), hash)
        hash = b'5d8a0231bb723e05696d6ef58a6da892'
        self._assert_hash(True, hash)
        self._assert_hash(1, hash)
        self._assert_hash(np.int(1), hash)
        self._assert_hash(np.int16(1), hash)
        self._assert_hash(np.int32(1), hash)
        self._assert_hash(1., hash)
        self._assert_hash(np.float(1), hash)
        self._assert_hash(np.float32(1), hash)
        self._assert_hash(np.float64(1), hash)

    def test_other(self):
        hash = b'083d6dde9ecc3001a8e931e81ba2ef5b'
        self._assert_hash(None, hash)

    def test_quantity(self):
        hash = b'815a7db5d24da2606b39c11790171428'
        self._assert_hash(ureg.Quantity([]), hash)
        hash = b'68da45d058243e374445ffa286ecec2f'
        self._assert_hash(ureg.Quantity(1), hash)
        hash = b'5207bcbf3207f96eb11a3d7b0a6220b9'
        self._assert_hash(ureg.Quantity([1, 2], 'mm'), hash)
        self._assert_hash(ureg.Quantity([1., 2], 'millimeter'), hash)
        self._assert_not_hash(ureg.Quantity([1.001, 2], 'millimeter'), hash)

    unittest.skip('not implemented yet')
    def test_uncertainties(self):
        pass

    def test_array(self):
        # Empty
        hash = b'425ed5b515f0a67c6316f0ec499f86bf'
        self._assert_hash(list(), hash)
        hash = b'bd2ce62877eb062903a7dca3b709c719'
        self._assert_hash(tuple(), hash)
        hash = b'c28c01d79239dc662f67b023595e6f91'
        self._assert_hash(set(), hash)
        hash = b'd1c2a7bc5a2517746e1788c5525d0dab'
        self._assert_hash(frozenset(), hash)
        hash = b'0336cebaf1e0ce25b6ac2f850e93dcc5'
        self._assert_hash(np.array([]), hash)

        # 0-dim
        hash = b'3ba6029fee48926bdc4542eaa1869081'
        self._assert_hash(np.array(5, dtype=np.int32), hash)
        self._assert_hash(np.array(5, dtype=np.float32), hash)

        # 1-dim
        hash = b'afe312eedebfd83274e9fc40c3125d22'
        self._assert_hash(list([1, 2]), hash)
        hash = b'afe312eedebfd83274e9fc40c3125d22'
        self._assert_hash(tuple([1, 2]), hash)
        hash = b'afe312eedebfd83274e9fc40c3125d22'
        self._assert_hash(set([1, 2]), hash)
        hash = b'afe312eedebfd83274e9fc40c3125d22'
        self._assert_hash(set([2, 1]), hash)
        hash = b'afe312eedebfd83274e9fc40c3125d22'
        self._assert_hash(frozenset([2, 1]), hash)
        hash = b'afe312eedebfd83274e9fc40c3125d22'
        self._assert_hash(np.array([1, 2]), hash)
        self._assert_hash(np.array([True, 2.]), hash)

    def test_dict(self):
        hash = b'73bc9620992249d67719f3170a967234'
        self._assert_hash({}, hash)
        hash = b'd34970312d6ed384fcc965f0a4d5d4f9'
        self._assert_hash({'a': 1}, hash)

    def _assert_hash(self, x, expected):
        self.assertEqual(hashing.calchash(x), expected)
        if not instance.isstring(x):
            self._assert_not_hash(str(x), expected)
            self._assert_not_hash(repr(x), expected)

    def _assert_not_hash(self, x, expected):
        self.assertNotEqual(hashing.calchash(x), expected)

    def _assert_hash_in(self, x, expected):
        hash = hashing.calchash(x)
        self.assertTrue(hash in expected, msg=hash)
        if not instance.isstring(x):
            self._assert_hash_not_in(str(x), expected)
            self._assert_hash_not_in(repr(x), expected)

    def _assert_hash_not_in(self, x, expected):
        self.assertFalse(hashing.calchash(x) in expected)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_hashing("test_string"))
    testSuite.addTest(test_hashing("test_number"))
    testSuite.addTest(test_hashing("test_other"))
    testSuite.addTest(test_hashing("test_array"))
    testSuite.addTest(test_hashing("test_quantity"))
    testSuite.addTest(test_hashing("test_uncertainties"))
    testSuite.addTest(test_hashing("test_dict"))
    testSuite.addTest(test_hashing("test_random"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

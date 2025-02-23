import unittest
import numpy as np
from .. import hashing
from ...patch.pint import ureg
from ...testutils.randomdata import factory
from .. import instance


class test_hashing(unittest.TestCase):
    def test_random(self):
        for _ in range(500):
            o = factory()
            a = o.data
            b = o.data
            _ = hashing.calchash(a)
            _ = hashing.calchash(b)
            self.assertTrue(hashing.hashequal(a, b))

    def test_string(self):
        # Empty string
        hash = b"7d52e4059c583de024b90a033ef4577e"
        self._assert_hash(b"", hash)
        self._assert_hash("", hash)
        # Ascii
        hash = b"4a77c5fa0df8eeb6f66ea3e95148e57e"
        self._assert_hash(b"abc", hash)
        self._assert_hash("abc", hash)
        # Extended ascii
        hash = b"9af537523552e3988d1eadc7dd4eb048"
        self._assert_hash(b"\xe423", hash)
        # UTF8
        hash = b"3a00931200bada8bb1a4c4c36185898f"
        self._assert_hash("\u0101bc", hash)

    def test_number(self):
        hash = b"cbeddc349fedb70a755c2b345e0f6a93"
        self._assert_hash(float("nan"), hash)
        hash = b"a713229c834747840d40cafa239c52ba"
        self._assert_hash(float("inf"), hash)
        hash = b"5d8a0231bb723e05696d6ef58a6da892"
        self._assert_hash(True, hash)
        self._assert_hash(1, hash)
        self._assert_hash(int(1), hash)
        self._assert_hash(np.int16(1), hash)
        self._assert_hash(np.int32(1), hash)
        self._assert_hash(1.0, hash)
        self._assert_hash(float(1), hash)
        self._assert_hash(np.float32(1), hash)
        self._assert_hash(np.float64(1), hash)

    def test_other(self):
        hash = b"083d6dde9ecc3001a8e931e81ba2ef5b"
        self._assert_hash(None, hash)

    def test_quantity(self):
        hash = b"6acad38f0f591c0434a820b6a5f3e4bf"
        self._assert_hash(ureg.Quantity([]), hash)
        hash = b"68da45d058243e374445ffa286ecec2f"
        self._assert_hash(ureg.Quantity(1), hash)
        hash = b"5207bcbf3207f96eb11a3d7b0a6220b9"
        self._assert_hash(ureg.Quantity([1, 2], "mm"), hash)
        self._assert_hash(ureg.Quantity([1.0, 2], "millimeter"), hash)
        self._assert_not_hash(ureg.Quantity([1.001, 2], "millimeter"), hash)

    unittest.skip("not implemented yet")

    def test_uncertainties(self):
        pass

    def test_array(self):
        # Empty
        hash = b"425ed5b515f0a67c6316f0ec499f86bf"
        self._assert_hash(list(), hash)
        hash = b"22b3afce62075c7012f8e5041adfee16"
        self._assert_hash(tuple(), hash)
        hash = b"f36acfd9d3040440cac6b32810abaeba"
        self._assert_hash(set(), hash)
        hash = b"3b9937ed5e5748f3a3081d71e017d706"
        self._assert_hash(frozenset(), hash)
        hash = b"e5906bb63f0e4c8da08aead693527d42"
        self._assert_hash(np.array([]), hash)

        # 0-dim
        hash = b"3ba6029fee48926bdc4542eaa1869081"
        self._assert_hash(np.array(5, dtype=np.int32), hash)
        self._assert_hash(np.array(5, dtype=np.float32), hash)

        # 1-dim
        hash = b"afe312eedebfd83274e9fc40c3125d22"
        self._assert_hash(list([1, 2]), hash)
        hash = b"afe312eedebfd83274e9fc40c3125d22"
        self._assert_hash(tuple([1, 2]), hash)
        hash = b"afe312eedebfd83274e9fc40c3125d22"
        self._assert_hash(set([1, 2]), hash)
        hash = b"afe312eedebfd83274e9fc40c3125d22"
        self._assert_hash(set([2, 1]), hash)
        hash = b"afe312eedebfd83274e9fc40c3125d22"
        self._assert_hash(frozenset([2, 1]), hash)
        hash = b"afe312eedebfd83274e9fc40c3125d22"
        self._assert_hash(np.array([1, 2]), hash)
        self._assert_hash(np.array([True, 2.0]), hash)

    def test_dict(self):
        hash = b"73bc9620992249d67719f3170a967234"
        self._assert_hash({}, hash)
        hash = b"d34970312d6ed384fcc965f0a4d5d4f9"
        self._assert_hash({"a": 1}, hash)

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

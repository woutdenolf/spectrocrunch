# -*- coding: utf-8 -*-

import os
import unittest
import itertools
import numpy as np
import h5py.h5t
from testfixtures import TempDirectory

from .. import h5fs
from .. import fs
from ..utils import TemporaryFilename
from ...testutils.subtest import TestCase


class test_h5fs(TestCase):
    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()

    def test_mode(self):
        # label   exist       not exist   read    write
        # r       -           error       -       error
        # r+      -           error       -       -
        # w       truncate    create      -       -
        # x       error       create      -       -
        # a       -           create      -       -
        #
        # Comparision with OS: w, x and a are always with '+' (i.e. read as well as write)

        with TemporaryFilename(self.dir.path, suffix=".txt") as filename:
            # Doesn't exist
            self.assertRaises((IOError, fs.Missing), self._check_r, filename, "123")
            self._check_w(filename, "123")
            # Exist
            self._check_r(filename, "456")
            self._check_content(filename, "123")

        with TemporaryFilename(self.dir.path, suffix=".txt") as filename:
            # Doesn't exist
            self.assertRaises((IOError, fs.Missing), self._check_rp, filename, "123")
            self._check_w(filename, "123")
            # Exist
            self._check_rp(filename, "456")
            self._check_content(filename, "123456")

        with TemporaryFilename(self.dir.path, suffix=".txt") as filename:
            # Doesn't exist
            self._check_w(filename, "123")
            self._check_content(filename, "123")
            # Exist
            self._check_w(filename, "456")
            self._check_content(filename, "456")

        with TemporaryFilename(self.dir.path, suffix=".txt") as filename:
            # Doesn't exist
            self._check_a(filename, "123")
            self._check_content(filename, "123")
            # Exist
            self._check_a(filename, "456")
            self._check_content(filename, "123456")

        with TemporaryFilename(self.dir.path, suffix=".txt") as filename:
            # Doesn't exist
            self._check_x(filename, "123")
            self._check_content(filename, "123")
            # Exist
            self.assertRaises(IOError, self._check_x, filename, "456")

    def _write(self, f, word):
        f.create_group(word)

    def _read(self, f, i=None):
        words = [k.split("/")[-1] for k in f]
        if i is not None:
            words = [words[-1]]
        return str("".join(words))

    def _check_w(self, filename, word):
        with h5fs.h5Device(filename, mode="w").open() as f:
            self._write(f, word)
            self.assertEqual(word, self._read(f))

    def _check_x(self, filename, word):
        with h5fs.h5Device(filename, mode="x").open() as f:
            self._write(f, word)
            self.assertEqual(word, self._read(f))

    def _check_a(self, filename, word):
        with h5fs.h5Device(filename, mode="a").open() as f:
            self._write(f, word)
            self.assertEqual(word, self._read(f, -1))

    def _check_r(self, filename, word):
        with h5fs.h5Device(filename, mode="r").open() as f:
            self.assertRaises(ValueError, self._write, f, word)
            self._read(f)

    def _check_rp(self, filename, word):
        with h5fs.h5Device(filename, mode="r+").open() as f:
            self._write(f, word)

    def _check_content(self, filename, word):
        filename = str(filename)
        if os.path.isfile(filename):
            with h5fs.h5Device(filename, mode="r").open() as f:
                b = self._read(f) == word
        else:
            b = None == word
        self.assertTrue(b)

    def test_path_splitting(self):
        cwd = self.dir.path
        func = lambda *args: os.path.join(cwd, *args)

        path = h5fs.Path(func("test.h5"))
        devsep = path.devsep

        self.assertEqual(path.device, func("test.h5"))
        self.assertEqual(path.path, "/")

        path = h5fs.Path(func("test.h5{}entry".format(devsep)))
        self.assertEqual(path.device, func("test.h5"))
        self.assertEqual(path.path, "/entry")

        path = h5fs.Path(func("test.h5{}/entry/abs{}def".format(devsep, devsep)))
        self.assertEqual(path.device, func("test.h5"))
        self.assertEqual(path.path, "/entry/abs{}def".format(devsep))

        path = h5fs.Path(func(".", "test.h5{}/entry/abs{}def".format(devsep, devsep)))
        self.assertEqual(path.device, os.path.abspath(func(".", "test.h5")))
        self.assertEqual(path.path, "/entry/abs{}def".format(devsep))

        path = h5fs.Path(
            func(".", "te{}st.h5{}/entry/abs{}def".format(devsep, devsep, devsep))
        )
        self.assertEqual(
            path.device, os.path.abspath(func(".", "te{}st.h5".format(devsep)))
        )
        self.assertEqual(path.path, "/entry/abs{}def".format(devsep))

        path = h5fs.Path(
            func(".", "te{}st{}/entry/abs{}def".format(devsep, devsep, devsep))
        )
        if path.sep == os.sep:
            self.assertEqual(
                path.device,
                os.path.abspath(
                    func(".", "te{}st{}/entry/abs{}def".format(devsep, devsep, devsep))
                ),
            )
            self.assertEqual(path.path, "/")
        else:
            self.assertEqual(
                path.device, os.path.abspath(func(".", "te{}st".format(devsep)))
            )
            self.assertEqual(path.path, "/entry/abs{}def".format(devsep))

        path = h5fs.Path(
            func(
                ".",
                "te{}st{}{}{}/entry/abs{}def".format(
                    devsep, devsep, devsep, devsep, devsep
                ),
            ),
            h5file=func(".", "te{}st".format(devsep)),
        )
        self.assertEqual(
            path.device, os.path.abspath(func(".", "te{}st".format(devsep)))
        )
        self.assertEqual(
            path.path, "/{}{}/entry/abs{}def".format(devsep, devsep, devsep)
        )

        path = h5fs.Path(
            func(
                ".",
                "te{}st{}{}{}/entry/abs{}def".format(
                    devsep, devsep, devsep, devsep, devsep
                ),
            ),
            h5file=func(".", "te{}st{}".format(devsep, devsep)),
        )
        self.assertEqual(
            path.device, os.path.abspath(func(".", "te{}st{}".format(devsep, devsep)))
        )
        self.assertEqual(path.path, "/{}/entry/abs{}def".format(devsep, devsep))

    def test_link_mixing(self):
        cwd = self.dir.path
        func = lambda *args: os.path.join(cwd, *args)

        roota = h5fs.Path(func("a.h5"))
        rootb = h5fs.Path(func("b.h5"))
        rootc = h5fs.Path(func("c.h5"))

        rootc["dir"].mkdir()
        fl = rootc["file.txt"].write(data=10)
        dest = rootc["dir"]["file.txt"].link(fl)
        a = roota["a"].link(rootb["b"])
        b = rootb["b"].link(dest)
        for f in [dest, a, b]:
            self.assertEqual(f.read(), 10)

        self.assertEqual(a.linkdest(), b)
        self.assertEqual(a.linkdest(follow=True), fl)
        self.assertEqual(b.linkdest(), dest)
        self.assertEqual(b.linkdest(follow=True), fl)
        self.assertEqual(dest.linkdest(), fl)

    def test_string(self):
        cwd = self.dir.path
        func = lambda *args: os.path.join(cwd, *args)
        for i, params in enumerate(itertools.product(*[(True, False)] * 3)):
            root = h5fs.Path(func("test_string{}.h5".format(i)))
            attribute, raiseExtended, useOpaqueDataType = params
            kwargs = {
                "attribute": attribute,
                "raiseExtended": raiseExtended,
                "useOpaqueDataType": useOpaqueDataType,
            }
            self._check_string(root, **kwargs)

    def _check_string(
        self, root, attribute=True, raiseExtended=True, useOpaqueDataType=True
    ):
        # Test following string literals
        sAsciiBytes = b"abc"
        sAsciiUnicode = "abc"
        sLatinBytes = b"\xe423"
        sLatinUnicode = "\xe423"  # not used
        sUTF8Unicode = "\u0101bc"
        sUTF8Bytes = b"\xc4\x81bc"
        sUTF8AsciiUnicode = "abc"
        sUTF8AsciiBytes = b"abc"
        # Expected conversion after HDF5 write/read
        strmap = {}
        strmap["ascii(scalar)"] = sAsciiBytes, sAsciiUnicode
        strmap["ext(scalar)"] = sLatinBytes, sLatinBytes
        strmap["unicode(scalar)"] = sUTF8Unicode, sUTF8Unicode
        strmap["unicode2(scalar)"] = sUTF8AsciiUnicode, sUTF8AsciiUnicode
        strmap["ascii(list)"] = (
            [sAsciiBytes, sAsciiBytes],
            [sAsciiUnicode, sAsciiUnicode],
        )
        strmap["ext(list)"] = [sLatinBytes, sLatinBytes], [sLatinBytes, sLatinBytes]
        strmap["unicode(list)"] = (
            [sUTF8Unicode, sUTF8Unicode],
            [sUTF8Unicode, sUTF8Unicode],
        )
        strmap["unicode2(list)"] = (
            [sUTF8AsciiUnicode, sUTF8AsciiUnicode],
            [sUTF8AsciiUnicode, sUTF8AsciiUnicode],
        )
        strmap["mixed(list)"] = (
            [sUTF8Unicode, sUTF8AsciiUnicode, sAsciiBytes, sLatinBytes],
            [sUTF8Bytes, sUTF8AsciiBytes, sAsciiBytes, sLatinBytes],
        )
        strmap["ascii(0d-array)"] = np.array(sAsciiBytes), sAsciiUnicode
        strmap["ext(0d-array)"] = np.array(sLatinBytes), sLatinBytes
        strmap["unicode(0d-array)"] = np.array(sUTF8Unicode), sUTF8Unicode
        strmap["unicode2(0d-array)"] = np.array(sUTF8AsciiUnicode), sUTF8AsciiUnicode
        strmap["ascii(1d-array)"] = (
            np.array([sAsciiBytes, sAsciiBytes]),
            [sAsciiUnicode, sAsciiUnicode],
        )
        strmap["ext(1d-array)"] = (
            np.array([sLatinBytes, sLatinBytes]),
            [sLatinBytes, sLatinBytes],
        )
        strmap["unicode(1d-array)"] = (
            np.array([sUTF8Unicode, sUTF8Unicode]),
            [sUTF8Unicode, sUTF8Unicode],
        )
        strmap["unicode2(1d-array)"] = (
            np.array([sUTF8AsciiUnicode, sUTF8AsciiUnicode]),
            [sUTF8AsciiUnicode, sUTF8AsciiUnicode],
        )
        strmap["mixed(1d-array)"] = (
            np.array([sUTF8Unicode, sUTF8AsciiUnicode, sAsciiBytes]),
            [sUTF8Unicode, sUTF8AsciiUnicode, sAsciiUnicode],
        )
        strmap["mixed2(1d-array)"] = (
            np.array([sUTF8AsciiUnicode, sAsciiBytes]),
            [sUTF8AsciiUnicode, sAsciiUnicode],
        )
        path = root.mkdir("test")

        prepare_kwargs = {
            "raiseExtended": raiseExtended,
            "useOpaqueDataType": useOpaqueDataType,
        }

        def write(name, value):
            if attribute:
                kwargs = {name: value}
                path.update_stats(prepare_kwargs=prepare_kwargs, **kwargs)
            else:
                path[name].mkfile(prepare_kwargs=prepare_kwargs, data=value)

        def read(name):
            if attribute:
                return path.get_stat(name)
            else:
                return path[name].read()

        # as long as vlen_opaque_dtype is not supported by h5py
        def remove00(s):
            return bytes(bytearray([b for b in bytearray(s) if b]))

        subtest_kwargs = {
            "attribute": attribute,
            "raiseExtended": raiseExtended,
            "useOpaqueDataType": useOpaqueDataType,
        }

        for name, (value, expectedValue) in strmap.items():
            subtest_kwargs["data"] = name
            with self.subTest(**subtest_kwargs):
                # Write/read
                decodingError = "ext" in name or name == "mixed(list)"
                expectOpaque = decodingError and useOpaqueDataType
                if raiseExtended and decodingError:
                    with self.assertRaises(UnicodeDecodeError):
                        write(name, value)
                    continue
                else:
                    write(name, value)
                value = read(name)

                # Check type and value
                if "list" in name or "1d-array" in name:
                    self.assertTrue(isinstance(value, np.ndarray))
                    if expectOpaque:
                        self.assertTrue(np.issubsctype(value.dtype, np.void))
                    value = value.tolist()  # also converts void to bytes
                    self.assertEqual(
                        list(map(type, value)), list(map(type, expectedValue)), msg=name
                    )
                    if expectOpaque:
                        value = list(map(remove00, value))
                    firstValue = value[0]
                else:
                    if expectOpaque:
                        value = remove00(value.tobytes())
                    firstValue = value
                msg = "{} {} instead of {}".format(
                    name, type(value), type(expectedValue)
                )
                self.assertEqual(type(value), type(expectedValue), msg=msg)
                self.assertEqual(value, expectedValue, msg=name)

                # Check HDF5 type id
                if not attribute:
                    with path[name].open() as dset:
                        typeId = dset.id.get_type()
                    if expectOpaque:
                        self.assertTrue(isinstance(typeId, h5py.h5t.TypeOpaqueID))
                    else:
                        self.assertTrue(isinstance(typeId, h5py.h5t.TypeStringID))
                        charSet = typeId.get_cset()
                        if isinstance(firstValue, bytes):
                            # This is why opaque types are used for extended ASCII strings
                            expectedCharSet = h5py.h5t.CSET_ASCII
                        else:
                            expectedCharSet = h5py.h5t.CSET_UTF8
                        msg = "{} type {} instead of {}".format(
                            name, charSet, expectedCharSet
                        )
                        self.assertEqual(charSet, expectedCharSet, msg=msg)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_h5fs("test_mode"))
    testSuite.addTest(test_h5fs("test_path_splitting"))
    testSuite.addTest(test_h5fs("test_link_mixing"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

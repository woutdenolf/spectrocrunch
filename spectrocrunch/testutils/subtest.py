# -*- coding: utf-8 -*-

import unittest
import sys
import itertools
from contextlib import contextmanager


class SkipTest(Exception):
    pass


class TestCase(unittest.TestCase):
    """Adds subTest for python 2"""

    hasSubTest = hasattr(unittest.TestCase, "subTest")

    def __init__(self, *args, **kwargs):
        self._inSubTest = False
        self._skipReason = ""
        super(TestCase, self).__init__(*args, **kwargs)

    @contextmanager
    def subTest(self, **kwargs):
        if self.hasSubTest:
            with super(TestCase, self).subTest(**kwargs):
                yield
            sys.stdout.write(".")
            sys.stdout.flush()
        else:
            self._inSubTest = True
            try:
                yield
            except SkipTest:
                pass
            finally:
                self._inSubTest = False
            sys.stdout.write(".")
            sys.stdout.flush()

    def skipTest(self, reason):
        if self._inSubTest:
            self._skipReason = reason
            raise SkipTest
        else:
            super(TestCase, self).skipTest(reason)

    @contextmanager
    def skipContext(self):
        self._skipReason = ""
        yield
        if self._skipReason:
            reason = self._skipReason
            self._skipReason = ""
            raise unittest.SkipTest(reason)

    def run_subtests(self, parameters, func):
        keys = list(parameters.keys())
        values = list(parameters.values())
        for ivalues in itertools.product(*values):
            kwargs = {key: value for key, value in zip(keys, ivalues)}
            with self.subTest(**kwargs):
                func(**kwargs)

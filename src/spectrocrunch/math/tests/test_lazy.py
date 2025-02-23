import unittest
import numpy as np
import operator
import random

from .. import lazy


class test_lazy(unittest.TestCase):
    def _gencase(self, ncases=100, logical=False):
        ops = [
            (operator.mul, 2),
            (operator.add, 2),
            (operator.sub, 2),
            (operator.pow, 2),
            (operator.truediv, 2),
            (operator.floordiv, 2),
            (operator.mod, 2),
            (abs, 1),
            (operator.pos, 1),
            (operator.neg, 1),
        ]

        if logical:
            ops += [
                (operator.lt, 2),
                (operator.le, 2),
                (operator.gt, 2),
                (operator.lt, 2),
                (operator.and_, 2),
                (operator.or_, 2),
                (operator.not_, 2),
                (operator.truth, 2),
            ]

        for i in range(ncases):
            n = random.randint(1, 50)
            yield [random.choice(ops) for i in range(n)]

    def _test_combine(self, logical=False):
        def randnr():
            return random.randint(0, 100) / 10.0 - 5

        def randbool():
            return random.choice([False, True])

        m = randnr()

        def randfunc():
            return lazy.Function(lambda x: m * x, "f(x0)")

        for ops in self._gencase(ncases=500, logical=logical):
            f = randfunc()
            x0 = randnr()

            x = f(x0)
            y = f

            for op, nargs in ops:
                if nargs == 1:
                    x = op(x)
                    y = op(y)
                elif nargs == 2:
                    if randbool():
                        a = randnr()
                        aeval = a
                    else:
                        a = randfunc()
                        aeval = a(x0)
                    right = randbool()

                    # Evaluate immediately
                    try:
                        if right:
                            x = op(x, aeval)
                        else:
                            x = op(aeval, x)
                    except Exception:
                        continue

                    # Evaluate Later
                    if right:
                        y = op(y, a)
                    else:
                        y = op(a, y)

            # Compare lazy evaluation with immediate result
            if np.isnan(x):
                self.assertTrue(np.isnan(y(x0)))
            else:
                self.assertEqual(x, y(x0))

    def test_arithm(self):
        self._test_combine(logical=False)

    def test_logical(self):
        self._test_combine(logical=True)

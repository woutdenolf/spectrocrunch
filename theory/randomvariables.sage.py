# -*- coding: utf-8 -*-

from sage.all_cmdline import *  # import sage library

_sage_const_2 = Integer(2)
_sage_const_1 = Integer(1)
_sage_const_0 = Integer(0)
_sage_const_1p5 = RealNumber("1.5")  #!/usr/bin/env sage

from sage.all import *
from sage.symbolic.integration.integral import definite_integral


def inverse(func, x):
    tmpvar = var("tmpvar", domain="real")
    return solve(x == func(tmpvar), tmpvar)[_sage_const_0].rhs()


def norm_pdf(x):
    return exp(-(x ** _sage_const_2) / _sage_const_2) / sqrt(_sage_const_2 * pi)


def norm_cdf(x):
    return (_sage_const_1 + erf(x / sqrt(_sage_const_2))) / _sage_const_2


def norm_ppf(x):
    return sqrt(_sage_const_2) * inverse(erf, _sage_const_2 * p - _sage_const_1)


def solvesinglevar(exp, var):
    seq = solve(exp, var)
    if len(seq) > _sage_const_1:
        raise RuntimeError("Expected one solution: {}".format(seq))
    eqn = seq[_sage_const_0]
    if eqn.lhs() != var:
        raise RuntimeError("Could not find solution to {}: {}".format(var, eqn))
    return eqn.rhs().simplify_full()


class rv(object):

    tmpvar = var("z", domain="real")

    def __init__(self, a=-infinity, b=infinity):
        self.a = a
        self.b = b

    def variable(self, name):
        v = var(name, domain="real")
        assume(self.a <= v <= self.b)
        return v

    def pvariable(self, name):
        v = var(name, domain="real")
        assume(_sage_const_0 <= v <= _sage_const_1)
        return v

    def pdf(self, x):
        raise NotImplementedError

    def cdf_from_pdf(self, x):
        return definite_integral(
            self.pdf(self.tmpvar), self.tmpvar, self.a, x
        ).simplify_full()

    def ppf_from_cdf(self, p):
        return solvesinglevar(p == self.cdf(self.tmpvar), self.tmpvar)

    def ppf(self, p):
        return self.ppf_from_cdf(p)

    def cdf(self, x):
        return self.cdf_from_pdf(x)

    def test(self):
        x = self.variable("x")
        p = self.pvariable("p")
        self._test(self.cdf(x), self.cdf_from_pdf(x))
        self._test(self.ppf(p), self.ppf_from_cdf(p))

    def _test(self, a, b):
        a = a.simplify_full()
        b = b.simplify_full()
        if not bool(a == b):
            raise RuntimeError("{} != {}".format(a, b))


class norm(rv):
    def pdf(self, x):
        return norm_pdf(x)

    def cdf(self, x):
        return norm_cdf(x)

    def ppf(self, p):
        return norm_ppf(p)


class truncnorm(rv):
    def __init__(self, a, b):
        super(truncnorm, self).__init__(a=a, b=b)
        self.cdfa = norm_cdf(a)
        self.cdfmina = norm_cdf(-a)
        self.cdfb = norm_cdf(b)
        self.cdfminb = norm_cdf(-b)
        if self.a > _sage_const_0:
            self.delta = -(self.cdfminb - self.cdfmina)
        else:
            self.delta = self.cdfb - self.cdfa

    def pdf(self, x):
        return norm_pdf(x) / self.delta

    def cdf(self, x):
        return (norm_cdf(x) - self.cdfa) / self.delta

    def ppf(self, p):
        if self.a > _sage_const_0:
            return -norm_ppf(p * self.cdfminb + self.cdfmina * (_sage_const_1 - p))
        else:
            return norm_ppf(p * self.cdfb + self.cdfa * (_sage_const_1 - p))


for o in [norm(), truncnorm(-_sage_const_1, _sage_const_1p5)]:
    o.test()

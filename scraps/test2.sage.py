# flake8: noqa

from sage.all_cmdline import *  # import sage library

_sage_const_1 = Integer(1)  #!/usr/bin/env sage

from sage.all import *


mu = function("mu")(z, E)

rho = function("rho")(z)

c = function("c")(z)

__tmp__ = var("x,y")
s = symbolic_expression(_sage_const_1 if y >= x else -_sage_const_1).function(x, y)

w = function("w")(i, z)

En = function("En")(i)

__tmp__ = var("x,y,a,b,E")
A = symbolic_expression(
    mu(a, E) * rho(a) * s(a, b) * (x - c(a))
    - mu(b, E) * rho(b) * s(a, b) * (y - c(b))
    - sum(mu(l, E) * rho(l), l, a + _sage_const_1, b - _sage_const_1)
).function(x, y, a, b, E)

__tmp__ = var("x,y,a,b,alpha,E")
T = symbolic_expression(exp(-A(x, y, a, b, E) / cos(alpha))).function(
    x, y, a, b, alpha, E
)

__tmp__ = var("a,b,i,j,z")
P = symbolic_expression(w(j, z) * rho(z) / cos(alpha)).function(a, b, i, j, z)

latex(P(a, b, i, j, z))

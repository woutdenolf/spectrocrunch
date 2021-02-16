# -*- coding: utf-8 -*-

from sage.all_cmdline import *  # import sage library

_sage_const_2 = Integer(2)
_sage_const_1 = Integer(1)
_sage_const_0 = Integer(0)
_sage_const_4 = Integer(4)  #!/usr/bin/env sage

from sage.all import *

# Number of layers
var("m")
assume(m, "integer")
assume(m > _sage_const_0)

# Number of interactions
var("n")
assume(n, "integer")
assume(n > _sage_const_0)

# Iterate over layers
var("l")
assume(l, "integer")
assume(l > _sage_const_0)
assume(l <= m)

# Absctract functions
mu = function("mu")
muF = function("muF")
rho = function("rho")
c = function("c")
__tmp__ = var("x,y")
s = symbolic_expression(_sage_const_1 if y >= x else -_sage_const_1).function(x, y)
w = function("w")
En = function("En")
d = function("d")

# Interaction cross-section
__tmp__ = var("i,E,theta,phi")
muint = symbolic_expression(muF(i) / (_sage_const_4 * pi)).function(i, E, theta, phi)

# Transmission
__tmp__ = var("x,y,a,b,E")
A = symbolic_expression(
    mu(a, E) * rho(a) * s(x, y) * (x - c(a))
    - mu(b, E) * rho(b) * s(x, y) * (y - c(b))
    - sum(mu(l, E) * rho(l) * d(l), l, a + _sage_const_1, b - _sage_const_1)
).function(x, y, a, b, E)

__tmp__ = var("x,y,a,b,alpha,E")
T = symbolic_expression(exp(-A(x, y, a, b, E) / cos(alpha))).function(
    x, y, a, b, alpha, E
)

# Interaction probablity
__tmp__ = var("j,E,z,alpha,theta,phi")
P = symbolic_expression(
    w(j, z) * rho(z) / cos(alpha) * muint(j, E, theta, phi)
).function(j, E, z, alpha, theta, phi)

# Number of interactions
n = _sage_const_1
J = [None] * (n + _sage_const_2)
J[_sage_const_0] = var("J_0")
E = list(var("E_%d" % i) for i in (ellipsis_iter(_sage_const_0, Ellipsis, n)))
z = list(var("z_%d" % i) for i in (ellipsis_iter(_sage_const_0, Ellipsis, n)))
z_0 = _sage_const_0
alpha = list(
    var("alpha_%d%d" % (i, i + _sage_const_1))
    for i in (ellipsis_iter(_sage_const_0, Ellipsis, n))
)

J[_sage_const_1] = (
    J[_sage_const_0]
    * T(x, y, alpha[_sage_const_0], E[_sage_const_0])
    * P(j, E[_sage_const_0], x, alpha[_sage_const_0], theta, phi)
).function(j, x, y, theta, phi)

print(J[_sage_const_1](j, z[_sage_const_0], z[_sage_const_1], theta, phi))

print(z_0)

# -*- coding: utf-8 -*-

from sage.all_cmdline import *  # import sage library

_sage_const_2 = Integer(2)
_sage_const_1 = Integer(1)
_sage_const_0 = Integer(0)  #!/usr/bin/env sage

from sage.all import *

# Phases
var("alpha,beta,delta,P", domain="real")

assume(P >= _sage_const_0, P <= _sage_const_1)

# Phasor sums
var("phix,phiy,phi", domain="complex")
assume(phi > _sage_const_0)

# Unpolarized:
ux = phix * exp(I * alpha)
uy = phiy * exp(I * alpha)
print(ux * conjugate(ux)).simplify_full()
print(ux * conjugate(uy)).simplify_full()
print(ux.real() * ux.real() + ux.imag() * ux.imag()).simplify_full().expand()

assert ux.real() * uy.real() + ux.imag() * uy.imag() == (ux * conjugate(uy)).real()
assert -ux.real() * uy.imag() + ux.imag() * uy.real() == (ux * conjugate(uy)).imag()

# Elliptical:
phix = phi
phiy = phi
cosbeta = cos(beta)
sinbeta = sin(beta)

# cosbeta = sqrt((1+P)/2)
# sinbeta = sqrt((1-P)/2)
sin2beta = _sage_const_2 * sinbeta * cosbeta
cos2beta = cosbeta ** _sage_const_2 - sinbeta ** _sage_const_2

ux0 = phix * cosbeta
uy0 = phiy * sinbeta
ux = ux0 * exp(I * alpha)
uy = uy0 * exp(I * delta) * exp(I * alpha)
# print (ux*conjugate(ux)).simplify_full()
# print (ux*conjugate(uy)).simplify_full()

uxr = ux.real()
uyr = uy.real()
A = phi.real()

AA = _sage_const_1 / (A * sin(delta) * cosbeta) ** _sage_const_2
CC = _sage_const_1 / (A * sin(delta) * sinbeta) ** _sage_const_2
BB = (
    -_sage_const_2
    * cos(delta)
    / (A ** _sage_const_2 * sin(delta) ** _sage_const_2 * sinbeta * cosbeta)
)

assert (
    AA * uxr ** _sage_const_2 + BB * uxr * uyr + CC * uyr ** _sage_const_2
    == _sage_const_1
)

DD = sqrt((AA - CC) ** _sage_const_2 + BB ** _sage_const_2)

semimajor = A * sqrt(
    (_sage_const_1 + sqrt(_sage_const_1 - (sin2beta * sin(delta)) ** _sage_const_2))
    / _sage_const_2
)
semiminor = A * sqrt(
    (_sage_const_1 - sqrt(_sage_const_1 - (sin2beta * sin(delta)) ** _sage_const_2))
    / _sage_const_2
)

tanchi = (CC - AA - DD) / BB
tan2chi = _sage_const_2 * tanchi / (_sage_const_1 - tanchi ** _sage_const_2)

assert tan2chi == sin2beta / cos2beta * cos(delta)

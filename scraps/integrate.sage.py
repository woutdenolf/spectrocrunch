# -*- coding: utf-8 -*-

from sage.all_cmdline import *  # import sage library

_sage_const_3 = Integer(3)
_sage_const_2 = Integer(2)
_sage_const_1 = Integer(1)
_sage_const_0 = Integer(0)
_sage_const_5 = Integer(5)
_sage_const_4 = Integer(4)
_sage_const_8 = Integer(8)  #!/usr/bin/env sage

from sage.all import *
from sage.symbolic.integration.integral import definite_integral

var("P")
var("theta")
var("phi")
var("R")
var("A")
var("B")
var("beta")
var("delta")

# Elastic scattering
Kunpol = (_sage_const_1 + (cos(theta)) ** _sage_const_2) / _sage_const_2
Kpol = Kunpol - (sin(theta)) ** _sage_const_2 / _sage_const_2 * (
    P * cos(_sage_const_2 * phi)
    + sqrt(_sage_const_1 - P**_sage_const_2) * cos(delta) * sin(_sage_const_2 * phi)
)
__tmp__ = var("theta,phi")
ElasticDiffPol2 = symbolic_expression(Kpol).function(theta, phi)
__tmp__ = var("theta")
ElasticDiffPol1 = symbolic_expression(
    definite_integral(
        ElasticDiffPol2(theta, phi), phi, _sage_const_0, _sage_const_2 * pi
    )
).function(theta)
__tmp__ = var("theta,phi")
ElasticDiffUnPol2 = symbolic_expression(Kunpol).function(theta, phi)
__tmp__ = var("theta")
ElasticDiffUnPol1 = symbolic_expression(
    definite_integral(
        ElasticDiffUnPol2(theta, phi), phi, _sage_const_0, _sage_const_2 * pi
    )
).function(theta)

assert ElasticDiffPol1(theta) == ElasticDiffUnPol1(theta)
assert (
    definite_integral(ElasticDiffUnPol1(theta) * sin(theta), theta, _sage_const_0, pi)
    == _sage_const_8 * pi / _sage_const_3
)

# Inelastic scattering
__tmp__ = var("theta,phi")
InelasticDiffPol2 = symbolic_expression(
    (_sage_const_1 / R + R + _sage_const_4 * Kpol - _sage_const_2)
    / (_sage_const_4 * R**_sage_const_2)
).function(theta, phi)
__tmp__ = var("theta")
InelasticDiffPol1 = symbolic_expression(
    definite_integral(
        InelasticDiffPol2(theta, phi), phi, _sage_const_0, _sage_const_2 * pi
    )
).function(theta)
__tmp__ = var("theta,phi")
InelasticDiffUnPol2 = symbolic_expression(
    (_sage_const_1 / R + R + _sage_const_4 * Kunpol - _sage_const_2)
    / (_sage_const_4 * R**_sage_const_2)
).function(theta, phi)
__tmp__ = var("theta")
InelasticDiffUnPol1 = symbolic_expression(
    definite_integral(
        InelasticDiffUnPol2(theta, phi), phi, _sage_const_0, _sage_const_2 * pi
    )
).function(theta)

assert InelasticDiffPol1(theta) == InelasticDiffUnPol1(theta)

print("\n" * _sage_const_5)
print("Azimuthally integrate differential Rayleigh cs:")
print(" Polarized source:", ElasticDiffPol1(theta))
print(" Unpolarized source:", ElasticDiffUnPol1(theta))
print("Differential Rayleigh cs with P=1 (synchrotron):")
print(
    ElasticDiffPol2(theta, phi)
    .substitute(P=_sage_const_1)
    .simplify_full()
    .simplify_trig()
)

print("")
print("Azimuthally integrate differential Compton cs:")
print(" Polarized source:", InelasticDiffPol1(theta))
print(" Unpolarized source:", InelasticDiffUnPol1(theta))
print("Differential Rayleigh cs with P=1 (synchrotron):")
print(InelasticDiffPol2(theta, phi).substitute(P=_sage_const_1).simplify_trig())

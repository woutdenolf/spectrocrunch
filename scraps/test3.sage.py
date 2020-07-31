# -*- coding: utf-8 -*-

from sage.all_cmdline import *  # import sage library

_sage_const_1 = Integer(1)
_sage_const_0 = Integer(0)  #!/usr/bin/env sage

from sage.all import *
from sage.symbolic.integration.integral import definite_integral

var("t,s,r,t", domain="real")
var("y,l", domain="integer")
var("ra,rb,rc", domain="real")
assume(ra > _sage_const_0)
assume(rb > _sage_const_0)
assume(rc > _sage_const_0)
assume(t > _sage_const_0)
assume(y >= _sage_const_0)
assume(l > _sage_const_0)

__tmp__ = var("t")
P0 = symbolic_expression(exp(-r * t)).function(t)
__tmp__ = var("t")
P = symbolic_expression(_sage_const_1 - P0(t)).function(t)
__tmp__ = var("t")
P1 = symbolic_expression(r * P0(t)).function(t)
__tmp__ = var("t")
P2 = symbolic_expression((P1(s) * P1(t - s)).integral(s, _sage_const_0, t)).function(t)
__tmp__ = var("t")
P3 = symbolic_expression((P2(s) * P1(t - s)).integral(s, _sage_const_0, t)).function(t)
__tmp__ = var("t")
P4 = symbolic_expression((P3(s) * P1(t - s)).integral(s, _sage_const_0, t)).function(t)


__tmp__ = var("t")
P1a = symbolic_expression(ra * exp(-ra * t)).function(t)
__tmp__ = var("t")
P1b = symbolic_expression(rb * exp(-rb * t)).function(t)
__tmp__ = var("t")
P1c = symbolic_expression(rb * exp(-rb * t)).function(t)
__tmp__ = var("t")
P2 = symbolic_expression((P1a(s) * P1b(t - s)).integral(s, _sage_const_0, t)).function(
    t
)
__tmp__ = var("t")
P3 = symbolic_expression((P2(s) * P1c(t - s)).integral(s, _sage_const_0, t)).function(t)
print P3

# g(y) = l^y/factorial(y)*exp(-l)
# print g.integral(y,0,oo)

#!/usr/bin/env sage

from sage.all import *
from sage.symbolic.integration.integral import definite_integral

var('t,s,r,t',domain='real')
var('y,l',domain='integer')
var('ra,rb,rc',domain='real')
assume(ra>0)
assume(rb>0)
assume(rc>0)
assume(t>0)
assume(y>=0)
assume(l>0)

P0(t) = exp(-r*t)
P(t) = 1-P0(t)
P1(t) = r*P0(t)
P2(t) = (P1(s)*P1(t-s)).integral(s,0,t)
P3(t) = (P2(s)*P1(t-s)).integral(s,0,t)
P4(t) = (P3(s)*P1(t-s)).integral(s,0,t)


P1a(t) = ra*exp(-ra*t)
P1b(t) = rb*exp(-rb*t)
P1c(t) = rb*exp(-rb*t)
P2(t) = (P1a(s)*P1b(t-s)).integral(s,0,t)
P3(t) = (P2(s)*P1c(t-s)).integral(s,0,t)
print P3

#g(y) = l^y/factorial(y)*exp(-l)
#print g.integral(y,0,oo)


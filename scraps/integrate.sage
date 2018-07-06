#!/usr/bin/env sage

from sage.all import *
from sage.symbolic.integration.integral import definite_integral

var('P')
var('theta')
var('phi')
var('R')
var('A')
var('B')
var('beta')
var('delta')

# Elastic scattering
Kunpol = (1+(cos(theta))**2)/2
Kpol = Kunpol - (sin(theta))**2/2*(P*cos(2*phi)+sqrt(1-P**2)*cos(delta)*sin(2*phi))
ElasticDiffPol2(theta,phi) = Kpol
ElasticDiffPol1(theta) = definite_integral(ElasticDiffPol2(theta,phi),phi,0,2*pi)
ElasticDiffUnPol2(theta,phi) = Kunpol
ElasticDiffUnPol1(theta) = definite_integral(ElasticDiffUnPol2(theta,phi),phi,0,2*pi)

assert(ElasticDiffPol1(theta)==ElasticDiffUnPol1(theta))
assert(definite_integral(ElasticDiffUnPol1(theta)*sin(theta),theta,0,pi)==8*pi/3)

# Inelastic scattering
InelasticDiffPol2(theta,phi) = (1/R + R + 4*Kpol - 2)/(4*R**2)
InelasticDiffPol1(theta) = definite_integral(InelasticDiffPol2(theta,phi),phi,0,2*pi)
InelasticDiffUnPol2(theta,phi) = (1/R + R + 4*Kunpol - 2)/(4*R**2)
InelasticDiffUnPol1(theta) = definite_integral(InelasticDiffUnPol2(theta,phi),phi,0,2*pi)

assert(InelasticDiffPol1(theta)==InelasticDiffUnPol1(theta))

print "\n"*5
print "Azimuthally integrate differential Rayleigh cs:"
print " Polarized source:",ElasticDiffPol1(theta)
print " Unpolarized source:",ElasticDiffUnPol1(theta)
print "Differential Rayleigh cs with P=1 (synchrotron):"
print ElasticDiffPol2(theta,phi).substitute(P=1).simplify_full().simplify_trig()

print ""
print "Azimuthally integrate differential Compton cs:"
print " Polarized source:",InelasticDiffPol1(theta)
print " Unpolarized source:",InelasticDiffUnPol1(theta)
print "Differential Rayleigh cs with P=1 (synchrotron):"
print InelasticDiffPol2(theta,phi).substitute(P=1).simplify_trig()




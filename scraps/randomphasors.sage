#!/usr/bin/env sage

from sage.all import *

# Phases
var('alpha,beta,delta,P',domain='real')

assume(P>=0,P<=1)

# Phasor sums
var('phix,phiy,phi',domain='complex')
assume(phi>0)

# Unpolarized:
ux = phix*exp(I*alpha)
uy = phiy*exp(I*alpha)
print (ux*conjugate(ux)).simplify_full()
print (ux*conjugate(uy)).simplify_full()
print (ux.real()*ux.real()+ux.imag()*ux.imag()).simplify_full().expand()

assert(ux.real()*uy.real()+ux.imag()*uy.imag()==(ux*conjugate(uy)).real())
assert(-ux.real()*uy.imag()+ux.imag()*uy.real()==(ux*conjugate(uy)).imag())

# Elliptical:
phix = phi
phiy = phi
cosbeta = cos(beta)
sinbeta = sin(beta)

#cosbeta = sqrt((1+P)/2)
#sinbeta = sqrt((1-P)/2)
sin2beta = 2*sinbeta*cosbeta
cos2beta = cosbeta**2-sinbeta**2

ux0 = phix*cosbeta
uy0 = phiy*sinbeta
ux = ux0*exp(I*alpha)
uy = uy0*exp(I*delta)*exp(I*alpha)
#print (ux*conjugate(ux)).simplify_full()
#print (ux*conjugate(uy)).simplify_full()

uxr = ux.real()
uyr = uy.real()
A = phi.real()

AA = 1/(A*sin(delta)*cosbeta)**2
CC = 1/(A*sin(delta)*sinbeta)**2
BB = -2*cos(delta)/(A**2*sin(delta)**2 *sinbeta*cosbeta)

assert( AA*uxr**2 + BB*uxr*uyr + CC*uyr**2 == 1) 

DD = sqrt((AA-CC)**2+BB**2)

semimajor = A*sqrt((1+sqrt(1-(sin2beta*sin(delta))**2))/2)
semiminor = A*sqrt((1-sqrt(1-(sin2beta*sin(delta))**2))/2)

tanchi = (CC-AA-DD)/BB
tan2chi = 2*tanchi/(1-tanchi**2)

assert( tan2chi==sin2beta/cos2beta*cos(delta))







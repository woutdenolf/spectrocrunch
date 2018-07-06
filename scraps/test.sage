#!/usr/bin/env sage

from sage.all import *

# Number of layers
var('m')
assume(m,'integer')
assume(m>0)

# Number of interactions
var('n')
assume(n,'integer')
assume(n>0)

# Iterate over layers
var('l')
assume(l,'integer')
assume(l>0)
assume(l<=m)

# Absctract functions
mu = function('mu')
muF = function('muF')
rho = function('rho')
c = function('c')
s(x,y) = 1 if y >= x else -1
w = function('w')
En = function('En')
d = function('d')

# Interaction cross-section
muint(i,E,theta,phi) = muF(i)/(4*pi)

# Transmission
A(x,y,a,b,E) = mu(a,E)*rho(a)*s(x,y)*(x-c(a)) - mu(b,E)*rho(b)*s(x,y)*(y-c(b)) - sum(mu(l,E)*rho(l)*d(l),l,a+1,b-1)

T(x,y,a,b,alpha,E) = exp(-A(x,y,a,b,E)/cos(alpha))

# Interaction probablity
P(j,E,z,alpha,theta,phi) = w(j,z)*rho(z)/cos(alpha)*muint(j,E,theta,phi)

# Number of interactions
n = 2
J = [None]*(n+1)
E = list(var('E_%d' % i) for i in (0..n))
z = list(var('z_%d' % i) for i in (0..n))
j = list(var('j_%d' % i) for i in (0..n))
alpha = list(var('alpha_%d%d' % (i,i+1)) for i in (0..n))

# Source
J[0] = var('J_0')
j[0] = None

# First interaction
J[1] = (J[0] * T(x,y,alpha[0],E[0]) * P(j[1],E[0],y,alpha[0],theta,phi)).function(x,y,theta,phi)

# Subsequent interactions
for k in range(2,n+1):
    J[k] = J[k-1](z[k-2],z[k-1],alpha[k-1],phi)*T(x,y,alpha[k-1],E[k-1])*P(j,E[k-1],z[k],alpha[k-1],theta,phi)


print 










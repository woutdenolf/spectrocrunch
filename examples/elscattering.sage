#!/usr/bin/env sage

from sage.all import *
var('theta,phi,S0,S1,S2,S3,psi,dolp,a,b,s,hdolp',domain='real')

S = vector([S0,S1,S2,S3])

def MuellerMatrixRotation(angle):
    return matrix(SR,4,4,[1,0,0,0,0,cos(2*angle),-sin(2*angle),0,0,sin(2*angle),cos(2*angle),0,0,0,0,1])
    
def JonesMatrixThomson(angle):
    return matrix(SR,2,2,[1,0,0,cos(angle)])

def jones_to_mueller(M):
    pauli44 = matrix(SR,4,4,[1,0,0,1,1,0,0,-1,0,1,1,0,0,1j,-1j,0])
    return (pauli44 * M.tensor_product(M.H) * pauli44.H)/2
 
M1 = MuellerMatrixRotation(phi)#.substitute(phi=pi/2).simplify()
M2 = jones_to_mueller(JonesMatrixThomson(theta)).substitute(theta = acos(sqrt(2*a-1))).simplify()
M = (M2*M1).simplify()
print(M)
S_scat = M*S

S_scat = S_scat.simplify()
S_scat = S_scat.substitute(S2=b*S1).simplify()
S_scat = S_scat.substitute(S1=s*(S0*dolp)/sqrt(1+b**2)).simplify()
print(S_scat[0])

#S_scat = S_scat.substitute(S2=b*S1).simplify()
#S_scat = S_scat.substitute(S1=hdolp*S0).simplify()
#print(S_scat[0])


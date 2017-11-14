#!/usr/bin/env sage

from sage.all import *
var('phasorsum,ux,uy,f',domain='complex')
var('alpha,beta,delta,P',domain='real')

var('theta0,varphi0,chi0,f0,ux0,uy0',domain='real')

var('theta1,varphi1,chi1,f1',domain='real')

var('theta2,varphi2,chi2,f2',domain='real')

assume(P>=0,P<=1)

ex0 = vector((1,0,0))
ey0 = vector((0,1,0))
ez0 = vector((0,0,1))
x0 = vector((cos(varphi0)*sin(theta0),sin(varphi0)*sin(theta0),cos(theta0)))
x1 = vector((cos(varphi1)*sin(theta1),sin(varphi1)*sin(theta1),cos(theta1)))
x2 = vector((cos(varphi2)*sin(theta2),sin(varphi2)*sin(theta2),cos(theta2)))

def Rot(v,urot,angle):
    return v*cos(angle)+urot.cross_product(v)*sin(angle)+urot.dot_product(v)*(1-cos(angle))*urot

def Re(v):
    return vector(x.real() for x in v)

def Im(v):
    return vector(x.imag() for x in v)

def Intensity(E):
    tmp = Re(E)
    return tmp.dot_product(tmp)

def Lx(angle):
    return Matrix([[0,0,0][0,cos(angle),-sin(angle)],[0,sin(angle),cos(angle)]])
    
def Ly(angle):
    return Matrix([[cos(angle),0,sin(angle)],[0,1,0],[-sin(angle),0,cos(angle)]])

def Lz(angle):
    return Matrix([[cos(angle),-sin(angle),0],[sin(angle),cos(angle),0],[0,0,1]])

def Simplify(x):
    try:
        return x.simplify_full().simplify_trig()
    except:
        return x

def SimplifyV(v):
    return vector(Simplify(x) for x in v)

def SimplifyM(m):
    return Matrix([SimplifyV(v) for v in m])

def Inner(v,w,e):
    return sum([v[i]*w[j]*e[i].dot_product(e[j]) for i in range(2) for j in range(2)])

def Scatter(xa,exa,eya,eza,Ea):
    uxa = Ea[0]
    uya = Ea[1]
    assert(Ea[2]==0)
    
    assert(exa.norm()==1)
    assert(eya.norm()==1)
    assert(eza.norm()==1)
    assert(xa.norm()==1)
    
    ezb = xa
    eyb = Simplify(eza.cross_product(ezb))
    exb = Simplify(eyb.cross_product(ezb))
    eyb = Simplify(eyb/eyb.norm())
    exb = Simplify(exb/exb.norm())
    
    assert(exb.dot_product(eyb)==0)
    assert(ezb.dot_product(eyb)==0)
    assert(ezb.dot_product(exb)==0)
    assert(Simplify(exb.norm())==1)
    assert(Simplify(eyb.norm())==1)
    assert(ezb.norm()==1)
    
    Eb = ((xa*exa)*xa-exa)*uxa + ((xa*eya)*xa-eya)*uya
    Eb = SimplifyV(Eb)
    print "\n   Scattered intensity:",Intensity(Eb)
    
    L = SimplifyM([exb,eyb,ezb])
    #C = Matrix([exb,eyb,ezb]).transpose()
    #L = SimplifyM(C.inverse())
    assert(SimplifyM(L.inverse())==SimplifyM(L.transpose()))
    Eb2 = SimplifyV(L*Eb)

    assert(Eb2[2]==0)
    assert(Simplify(Intensity(Eb))==Simplify(Intensity(Eb2)))
    
    return exb,eyb,ezb,Eb2

#cosbeta = sqrt((1+P)/2)
#sinbeta = sqrt((1-P)/2)
#cosbeta = cos(beta)
#sinbeta = sin(beta)
#ux = phasorsum*cosbeta*exp(I*alpha)
#uy = phasorsum*sinbeta*exp(I*delta)*exp(I*alpha)

Sx = ux*f
Sy = uy*f
print((Sx.real())^2).expand()
print(Sx.real()*Sy.real()).expand()

#E0 = ux0*ex0 + uy0*ey0
#print "\nSource intensity:",Intensity(E0)

#ex1,ey1,ez1,E1 = Scatter(x0,ex0,ey0,ez0,E0)
#print "\nScattered intensity:",Simplify(Intensity(E1))

#ex2,ey2,ez2,E2 = Scatter(x1,ex1,ey1,ez1,E1)
#print "\nScattered intensity:",Simplify(Intensity(E2))




#!/usr/bin/env sage

from sage.all import *
from sage.symbolic.integration.integral import definite_integral

def MuellerMatrixRotation(angle):
        return matrix(SR,4,4,[1,0,0,0,0,cos(2*angle),-sin(2*angle),0,0,sin(2*angle),cos(2*angle),0,0,0,0,1])
        
def JonesMatrixThomson(angle):
    return matrix(SR,2,2,[1,0,0,cos(angle)])

pauli44 = matrix(SR,4,4,[1,0,0,1,1,0,0,-1,0,1,1,0,0,1j,-1j,0])
pauli44i = pauli44.H/2

def jones_to_mueller(M):
    return pauli44 * M.tensor_product(M.H) * pauli44i

def mueller_to_jones(M):
    return pauli44i * M * pauli44 # decompose in M.tensor_product(M.H) not generally possible

def thomson():
    var('theta,phi,beta,S0,S1,S2,S3,a',domain='real')
    
    S = vector([S0,S1,S2,S3])

    M1 = MuellerMatrixRotation(beta)
    M2 = jones_to_mueller(JonesMatrixThomson(theta))
    #print mueller_to_jones(M2)
    
    M2 = M2.substitute(theta = acos(sqrt(2*a-1))).simplify()
    M = (M2*M1).simplify()
    print "\nMueller matrix:"
    print M
    S_scat = (M*S).simplify()

    I_scat = S_scat[0].substitute(beta=pi/2-phi).simplify()
    print "\nScattered intensity:"
    print I_scat

    It_scat = I_scat.substitute(a=(1+cos(theta)**2)/2)
    
    print "\nScattered intensity (unpol):"
    print It_scat.substitute(S1=0).substitute(S2=0).simplify_full()

    print "\nScattered intensity (horizontal linear pol.):"
    print It_scat.substitute(S1=S0).substitute(S2=0).simplify_full()

    It_scat = definite_integral(It_scat*sin(theta),theta,0,pi)
    It_scat = definite_integral(It_scat,phi,0,2*pi)
    print "\nThomson cross-section (units of r_e^2):"
    print It_scat/S0

def compton():
    var('theta,phi,beta,S0,S1,S2,S3,a,c,costheta,s,E,Esc',domain='real')
    assume(E>0)
    assume(Esc>0)  
    
    S = vector([S0,S1,S2,S3])
    
    M1 = MuellerMatrixRotation(beta)
    M2 = matrix(SR,4,4,[a+c,1-a,0,0,1-a,a,0,0,0,0,cos(theta),0,0,0,0,cos(theta)*(1+c)])*s
    #print mueller_to_jones(M2)
    M2 = M2.substitute(theta = acos(sqrt(2*a-1))).simplify()
    
    M = (M2*M1).simplify()
    print "\nMueller matrix:"
    print M
    S_scat = (M*S).simplify()
    
    I_scat = S_scat[0].substitute(beta=pi/2-phi).simplify()
    print "\nScattered intensity:"
    print I_scat
    
    # Energy in units of m_e*c^2
    It_scat = I_scat.substitute(a=(1+cos(theta)**2)/2)
    It_scat = It_scat.substitute(s=Esc^2/E^2)
    It_scat = It_scat.substitute(c=(E-Esc)/2*(1-cos(theta)))
    It_scat = It_scat.simplify_full()
    
    print "\nScattered intensity (unpol):"
    print It_scat.substitute(S1=0).substitute(S2=0).simplify_full()
    
    It_scatu1 = It_scat.substitute(S1=0).substitute(S2=0)
    It_scatu2 = Esc**2/E**2*S0*(E/Esc+Esc/E-sin(theta)^2)/2
    It_scatu1 = It_scatu1.substitute(Esc=E/(1+E*(1-cos(theta)))).simplify_full().simplify_trig()
    It_scatu2 = It_scatu2.substitute(Esc=E/(1+E*(1-cos(theta)))).simplify_full().simplify_trig()
    print bool(It_scatu1==It_scatu2)
    
    It_scatu1 = It_scat.substitute(S1=S0).substitute(S2=0)
    It_scatu2 = Esc**2/E**2*S0*(E/Esc+Esc/E-2*sin(theta)^2*cos(phi)^2)/2
    It_scatu1 = It_scatu1.substitute(Esc=E/(1+E*(1-cos(theta)))).simplify_full().simplify_trig()
    It_scatu2 = It_scatu2.substitute(Esc=E/(1+E*(1-cos(theta)))).simplify_full().simplify_trig()
    print bool(It_scatu1==It_scatu2)

    print "\nScattered intensity (horizontal linear pol.):"
    print It_scat.substitute(S1=S0).substitute(S2=0).simplify_full()

    It_scat = It_scat.substitute(Esc=E/(1+E*(1-cos(theta))))
    It_scat = It_scat.simplify_full()
    It_scat = definite_integral(It_scat*sin(theta),theta,0,pi)
    It_scat = definite_integral(It_scat,phi,0,2*pi)
    print "\nKlein-Nishina cross-section (units of r_e^2):"
    print (It_scat/S0).simplify_full()
       
thomson()
compton()



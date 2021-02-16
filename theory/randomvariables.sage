#!/usr/bin/env sage

from sage.all import *
from sage.symbolic.integration.integral import definite_integral

def inverse(func,x):
    tmpvar = var('tmpvar',domain='real')
    return solve( x == func(tmpvar), tmpvar)[0].rhs()

def norm_pdf(x):
    return exp(-x**2/2)/sqrt(2*pi)

def norm_cdf(x):
    return (1+erf(x/sqrt(2)))/2

def norm_ppf(x):
    return sqrt(2)*inverse(erf,2*p - 1)

def solvesinglevar(exp,var):
    seq = solve(exp,var)
    if len(seq)>1:
        raise RuntimeError("Expected one solution: {}".format(seq))
    eqn = seq[0]
    if eqn.lhs()!=var:
        raise RuntimeError("Could not find solution to {}: {}".format(var,eqn))
    return eqn.rhs().simplify_full()

class rv(object):

    tmpvar = var('z',domain='real')

    def __init__(self,a=-infinity,b=infinity):
        self.a = a
        self.b = b

    def variable(self,name):
        v = var(name,domain='real')
        assume(self.a<=v<=self.b)
        return v
    
    def pvariable(self,name):
        v = var(name,domain='real')
        assume(0<=v<=1)
        return v
        
    def pdf(self,x):
        raise NotImplementedError

    def cdf_from_pdf(self,x):
        return definite_integral(self.pdf(self.tmpvar),self.tmpvar,self.a,x).simplify_full()

    def ppf_from_cdf(self,p):
        return solvesinglevar(p==self.cdf(self.tmpvar),self.tmpvar)
        
    def ppf(self,p):
        return self.ppf_from_cdf(p)
        
    def cdf(self,x):
        return self.cdf_from_pdf(x)

    def test(self):
        x = self.variable('x')
        p = self.pvariable('p')
        self._test(self.cdf(x),self.cdf_from_pdf(x))
        self._test(self.ppf(p),self.ppf_from_cdf(p))

    def _test(self,a,b):
        a = a.simplify_full()
        b = b.simplify_full()
        if not bool(a==b):
            raise RuntimeError("{} != {}".format(a,b))

class norm(rv):

    def pdf(self,x):
        return norm_pdf(x)

    def cdf(self,x):
        return norm_cdf(x)

    def ppf(self,p):
        return norm_ppf(p)

class truncnorm(rv):

    def __init__(self,a,b):
        super(truncnorm,self).__init__(a=a,b=b)
        self.cdfa = norm_cdf(a)
        self.cdfmina = norm_cdf(-a)
        self.cdfb = norm_cdf(b)
        self.cdfminb = norm_cdf(-b)
        if self.a>0:
            self.delta = -(self.cdfminb - self.cdfmina)
        else:
            self.delta = self.cdfb - self.cdfa
            
    def pdf(self,x):
        return norm_pdf(x)/self.delta

    def cdf(self,x):
        return (norm_cdf(x) - self.cdfa)/self.delta

    def ppf(self,p):
        if self.a>0:
            return -norm_ppf(p*self.cdfminb + self.cdfmina*(1-p))
        else:
            return norm_ppf(p*self.cdfb + self.cdfa*(1-p))

for o in [norm(),truncnorm(-1,1.5)]:
    o.test()


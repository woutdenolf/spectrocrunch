# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import numpy as np
import scipy.optimize
import warnings

def gaussian(x,x0,sx,A):
    return A/(np.sqrt(2*np.pi)*sx)*np.exp(-(x-x0)**2/(2*sx**2))

def guess_gaussian(x,data):
    x0i = np.argmax(data)
    x0 = x[x0i]
    sx = np.sqrt(abs((x-x0)**2*data).sum()/data.sum())
    A = data[x0]*np.sqrt(2*np.pi)*sx
    return np.array([x0,sx,A],dtype=np.float32)

def fitgaussian(x,data):
    return leastsq(x,data,guessfunc=guess_gaussian,func=gaussian)

def leastsq(x,data,guessfunc=None,fitfunc=None):
    guess = guessfunc(x,data)
    errorfunc = lambda p,x,data: np.ravel(fitfunc(x,*p)-data)
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p, success = scipy.optimize.leastsq(errorfunc, guess, args=(x,data))
        success = success>0 and success<5

    return p, success

def xyremovenan(x,y):
    b = np.logical_and(~np.isnan(x),~np.isnan(y))
    x[b],y[b]

def cor_from_cov(cov):
    D = np.diag(1/np.sqrt(np.diag(cov)))
    return D.dot(cov.dot(D))
    
def lstsq_cov(A,vare=None,cove=None):
    # A.x = b + e
    # E(e) = 0
    #
    # x = A^(-1).(b+e)
    # COVx = A^(-1).COVe.A^(-T)
    # COVx = (A^T.COVe^(-1).A)^(-1)
    
    #try:
    if cove is None:
        return np.linalg.inv(np.dot(A.T, A/vare.reshape((vare.size,1))))
    else:
        iA = np.linalg.inv(A)
        return invA.dot(COVe.dot(invA.T))
    #except np.linalg.linalg.LinAlgError:
    #    _,n = A.shape
    #    return np.ones((n,n))
    
def lstsq_cov_est(A,b,x):
    # A.x = b + e
    # E(e) = 0
    #
    # x = A^(-1).(b+e)
    # COVx = A^(-1).COVe.A^(-T)
    # COVx = (A^T.COVe^(-1).A)^(-1)
    #
    # COVe diagonal (with VARe on the diagonal)
    # VAR(e) ≃ VAR(A.x-b)
    
    vare = np.var(np.dot(A, x)-b, ddof=x.size)
    return np.linalg.inv(np.dot(A.T, A)) * vare

def lstsq_std(A,b=None,x=None,vare=None):
    """Estimated error of solution to linear system
    
    .. math::
    
        A.x = b + e
        
        E(e) = 0
    
    Args:
        A(array): (m x n)
        b(array): (m)
        x(array): (n)
        vare(array): variance on e (m)
        
    Returns:
        stdx(array): errors (n)
    """
    if vare is None:
        return np.sqrt(np.diag(lstsq_cov_est(A,b,x)))
    else:
        return np.sqrt(np.diag(lstsq_cov(A,vare)))

def lstsq_std_indep(A,b=None,x=None,vare=None):
    """Estimated error of solution to linear system
    
    .. math::
    
        A.x = b + e
        
        E(e) = 0
    
    Assume we know x are independent random variables then
    there variances are found by solving another linear system
    
    .. math::
        (A*A).VARx = VARe
    
    Args:
        A(array): (m x n)
        b(array): (m)
        x(array): (n)
        vare(array): variance on e (m)
        
    Returns:
        stdx(array): sqrt(VARX) (n)
    """
    if vare is None:
        vare = np.var(np.dot(A, x)-b, ddof=x.size)
        vare = np.full(vare,x.size)
    return np.sqrt(lstsq(A*A,vare))

def lstsq(A,b,errors=False,vare=None):
    """Solve the following linear system
    
    .. math::
    
        A.x = b + e
        
        E(e) = 0
    
    Args:
        A(array): (m x n)
        b(array): (m)
        errors(Optional(bool)): return solution with estimated error
        vare(array): variance on e (m)
        
    Returns:
        x(array): solution (n)
        stdx(array): optional errors (n)
    """
    x = np.linalg.lstsq(A, b, rcond=-1)[0]
    if errors:
        return x,lstsq_std(A,b=b,x=x,vare=vare)
    else:
        return x

def lstsq_nonnegative(A,b,errors=False,vare=None):
    """Solve the following linear system
    
    .. math::
    
        A.x = b + e \quad x>=0
        
        E(e) = 0
    
    Args:
        A(array): (m x n)
        b(array): (m)
        errors(Optional(bool)): return solution with estimated error
        vare(array): variance on e (m)
        
    Returns:
        x(array): solution (n)
        stdx(array): optional errors (n)
    """
    x = scipy.optimize.nnls(A, b)[0]
    if errors:
        return x,lstsq_std(A,b=b,x=x,vare=vare)
    else:
        return x
        
def lstsq_bound(A,b,lb,ub,errors=False,vare=None):
    """Solve the following linear system
    
    .. math::
    
        A.x = b + e \quad lb<=x<=ub
        
        E(e) = 0
    
    Args:
        A(array): (m x n)
        b(array): (m)
        lb(num): lower bound
        ub(num): upper bound
        errors(Optional(bool)): return solution with estimated error
        vare(array): variance on e (m)
        
    Returns:
        x(array): solution (n)
        stdx(array): optional errors (n)
    """
    x = scipy.optimize.lsq_linear(A, b, bounds=(lb,ub)).x
    if errors:
        return x,lstsq_std(A,b=b,x=x,vare=vare)
    else:
        return x
        
def linfit(x,y,errors=False,vare=None):
    """Linear fit
    
    .. math::
    
        y = m.x + b + e
        
        E(e) = 0
    
    Args:
        x(array): (m)
        y(array): (m)
        errors(Optional(bool)): return solution with estimated error
        vare(array): variance on e (m)
        
    Returns:
        sol(array): solution (m,b)
        stdsol(array): optional errors (2-tuple)
    """
    A = np.vstack([x, np.ones(len(x))]).T
    return lstsq(A,y,errors=errors,vare=vare) # slope,intercept

def linfit2(x,y,errors=False,vare=None):
    if vare is not None:
        raise NotImplementedError("Use linfit instead")
    n = len(x)
    Sxy = (x*y).sum()
    Sxx = (x*x).sum()
    Sx = x.sum()
    Sy = y.sum()
    denom = float(n*Sxx-Sx*Sx)
    mnum = n*Sxy-Sx*Sy
    bnum = Sxx*Sy-Sx*Sxy
    m = mnum / denom
    b = bnum / denom
    
    if errors:
        Syy = (y*y).sum()
        num = n*Syy-Sy*Sy-m*mnum
        mstd = np.sqrt(    num / ((n-2.)*denom)  )
        bstd = np.sqrt(num*Sxx / (n*(n-2.)*denom))
        return [m,b],[mstd,bstd]
    else:
        return [m,b]

def nanlinfit(x,y,errors=False,vare=None):
    x,y = xyremovenan(x,y)
    return linfit(x,y,errors=errors,vare=vare)

def nanlinfit2(x,y,errors=False,vare=None):
    x,y = xyremovenan(x,y)
    return linfit2(x,y,errors=errors,vare=vare)

def linfit_zerointercept(x,y,errors=False,vare=None):
    """Linear fit with zero intercept
    
    .. math::
    
        y = m.x + e
        
        E(e) = 0
    
    Args:
        x(array): (m)
        y(array): (m)
        errors(Optional(bool)): return solution with estimated error
        vare(array): variance on e (m)
        
    Returns:
        m(array): solution
        stdm(array): optional error
    """
    A = np.vstack([x]).T
    if errors:
        m,mstd = lstsq(A,y,errors=True,vare=vare)
        return m[0],mstd[0]
    else:
        return lstsq(A,y)[0]
    
def linfit_zerointercept2(x,y,errors=False,vare=None):
    Sxy = (x*y).sum()
    Sxx = float((x*x).sum())
    m = Sxy/Sxx
    if errors:
        n = len(x)
        Syy = (y*y).sum()
        mstd = np.sqrt( (Syy+m*m*Sxx-2*m*Sxy) / ((n-1.)*Sxx) ) # Not sure
        return m,mstd
    return m


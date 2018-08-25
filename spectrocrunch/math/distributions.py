# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
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

from scipy.stats import rv_continuous
import scipy.special as special
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import _continuous_distns as rv_helper

def plothistogram(values,edges=None,**kwargs):
    """Usage:
        plothistogram(samples,density=1)
        plothistogram(*np.histogram(samples,density=1))
       
       Args:
        values(array): sample or histogram (when edges is set)
        edges(Optional(array)): histrogram bin edges
    """
    if edges is None:
        plt.hist(values,align='mid',alpha=0.7,**kwargs)
    else:
        wbin = np.diff(edges)
        xleft = edges[:-1]
        plt.bar(xleft,height=values,width=wbin,align='edge',**kwargs)
        xcen = xleft + wbin*0.5
        plt.plot(xcen,values,'o') 
        return xcen


class limitednorm_gen(rv_continuous):
    """Normal distribution within finite limits [-k,k]
    
       pdf(x) = ptf_norm(x) + cdf_norm(-k)
    """

    def _pdf(self,y):
        # y = (x-loc)/scale

        invalid = (y<-self.k) | (y>self.k)
        if isinstance(y,np.ndarray):
            ret = np.exp(-y**2 / 2.)/np.sqrt(2.0 * np.pi) + self.offset
            ret[invalid] = 0
            return ret
        elif invalid:
            return 0*y
        else:
            return np.exp(-y**2 / 2.)/np.sqrt(2.0 * np.pi) + self.offset

    def _cdf(self,y):
        # y = (x-loc)/scale

        if isinstance(y,np.ndarray):
            ret = special.ndtr(y) + self.offset
            ret[y<-self.k] = 0
            ret[y>self.k] = 1
            return ret
        else:
            if y<=-self.k:
                return 0
            elif y>=self.k:
                return 1
            else:
                return special.ndtr(y) + self.offset

    def _ppf(self,p):
        return special.ndtri(p-self.offset)


class truncnorm_gen(rv_continuous):
    """Normal distribution within finite limits [a,b]
    """
    
    def _argcheck(self, a, b):
        self.a = a
        self.b = b
        self._cdfb = rv_helper._norm_cdf(b)
        self._cdfa = rv_helper._norm_cdf(a)
        self._cdfminb = rv_helper._norm_cdf(-b)
        self._cdfmina = rv_helper._norm_cdf(-a)
        self._delta = np.where(self.a > 0,
                               -(self._cdfminb - self._cdfmina),
                               self._cdfb - self._cdfa)
        self._logdelta = np.log(self._delta)
        return a != b

    def _pdf(self, x, a, b):
        return rv_helper._norm_pdf(x) / self._delta

    def _logpdf(self, x, a, b):
        return rv_helper._norm_logpdf(x) - self._logdelta

    def _cdf(self, x, a, b):
        return (rv_helper._norm_cdf(x) - self._cdfa) / self._delta

    def _ppf(self, q, a, b):
        return np.where(self.a > 0,
                       -rv_helper._norm_ppf(q*self._cdfminb + self._cdfmina*(1.0-q)),
                        rv_helper._norm_ppf(q*self._cdfb + self._cdfa*(1.0-q)))

    def _stats(self, a, b):
        nA, nB = self._cdfa, self._cdfb
        d = nB - nA
        pA, pB = rv_helper._norm_pdf(a), rv_helper._norm_pdf(b)
        mu = (pA - pB) / d   # correction sign
        mu2 = 1 + (a*pA - b*pB) / d - mu*mu
        return mu, mu2, None, None


def limitednorm(k,**kwargs):
    k = np.abs(k)
    return truncnorm_gen(name='limitednorm',**kwargs)(a=-k,b=k)
    #return scipy.stats.truncnorm(a=-k,b=k)
    
    
class holenorm_gen(rv_continuous):
    """Flipped normal distribution within finite limits [a,b]
    """
    
    def _argcheck(self, a, b):
        self.a = a
        self.b = b
        self._cdfb = rv_helper._norm_cdf(b)
        self._cdfa = rv_helper._norm_cdf(a)
        self._cdfminb = rv_helper._norm_cdf(-b)
        self._cdfmina = rv_helper._norm_cdf(-a)
        self._delta = np.where(self.a > 0,
                               -(self._cdfminb - self._cdfmina),
                               self._cdfb - self._cdfa)
        self._logdelta = np.log(self._delta)
        return a != b

    def _pdf(self, x, a, b):
        return rv_helper._norm_pdf(x) / self._delta

    def _logpdf(self, x, a, b):
        return rv_helper._norm_logpdf(x) - self._logdelta

    def _cdf(self, x, a, b):
        return (rv_helper._norm_cdf(x) - self._cdfa) / self._delta

    def _ppf(self, q, a, b):
        return np.where(self.a > 0,
                       -rv_helper._norm_ppf(q*self._cdfminb + self._cdfmina*(1.0-q)),
                        rv_helper._norm_ppf(q*self._cdfb + self._cdfa*(1.0-q)))

    def _stats(self, a, b):
        nA, nB = self._cdfa, self._cdfb
        d = nB - nA
        pA, pB = rv_helper._norm_pdf(a), rv_helper._norm_pdf(b)
        mu = (pA - pB) / d   # correction sign
        mu2 = 1 + (a*pA - b*pB) / d - mu*mu
        return mu, mu2, None, None


def holenorm(k,**kwargs):
    k = np.abs(k)
    return holenorm_gen(name='holenorm',**kwargs)(a=-k,b=k)
    
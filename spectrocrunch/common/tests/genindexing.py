# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
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

import random

import itertools

import numpy as np

def genindexing(dim,advanced=False):

    if dim>1:
        i = random.randint(1,dim-1)
        j = random.randint(-dim,-1)
        k = random.randint(1,max(1,dim//2))
    else:
        i = 0
        j = 0
        k = 1

    ret = [0,i,j,\
           slice(i,j,k),slice(None,j,k),slice(i,None,k),slice(i,j,None),\
           slice(None,None,k),slice(i,None,None),slice(None,j,None),\
           slice(None),\
           slice(0)]

    # Advanced indexing
    if advanced:
        ret += [[random.randint(-dim,dim-1),random.randint(-dim,dim-1)],\
                [bool(random.getrandbits(1))]*dim,[True]*dim,[False]*dim,\
                Ellipsis,np.newaxis]

    return ret

def valid(index):
    # No more than one '...'
    if sum([i is Ellipsis for i in index])>1:
        return False

    # More than one list: they must have the same length
    b = [i for i in index if isinstance(i,list)]
    if len(b)>1:
        b = [len(i) for i in b]
        if b.count(b[0]) != len(b):
            return False

    return True

def genindexingn(dim,advanced=False):
    lst = [genindexing(n,advanced=advanced) for n in dim]
    lst = list(itertools.product(*lst))
    return [t for t in lst if valid(t)]

        

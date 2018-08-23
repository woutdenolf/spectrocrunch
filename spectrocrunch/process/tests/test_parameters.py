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

import unittest

from .. import parameters
from ...common import hashing
import numpy as np
import random
import collections

class test_parameters(unittest.TestCase):

    def hashequal(self,a,b,**kwargs):
        self.assertEqual(hashing.calchash(a,**kwargs),hashing.calchash(b,**kwargs))

    def hashnotequal(self,a,b,**kwargs):
        self.assertNotEqual(hashing.calchash(a,**kwargs),hashing.calchash(b,**kwargs))

    def test_hash(self):
        f1 = 1.3
        f2 = 2.4
        
        # Numbers:
        self.hashnotequal(np.float32(f1),np.float64(f1))
        self.hashequal(float("1.23e-4"),1.23e-4)

        # Dictionaries
        self.hashequal({'a':f1,'b':f2},{'a':f1,'b':f2})
        self.hashnotequal(collections.OrderedDict(zip(['a','b'],[f1,f2])),collections.OrderedDict(zip(['b','a'],[f2,f1])))
        
        # Numpy arrays
        a = np.arange(10*20*30).reshape((10,20,30))
        b = np.arange(10*20*30).reshape((10,20,30))
        self.hashequal((1,a),(1,b))
        self.hashequal(a,b.tolist())
        self.hashnotequal((1,a),(1,b.tolist()),numpylarge=True)

        # Nesting of dictionaries and sequences
        ran1 = np.linspace(0,f1,2)
        ran2 = np.linspace(0,f2,3)
        seq = [lambda x:set(sorted(x, key=lambda k: random.random())),
               lambda x:frozenset(sorted(x, key=lambda k: random.random())),
               tuple,
               np.array,
               lambda x:x]
        for seq1 in seq:
            for seq2 in seq:
                d1 = dict(zip(['a','b'],[seq1(ran1),seq2(ran2)]))
                d2 = dict(zip(['b','a'],[seq1(ran2),seq2(ran1)]))
                d3 = dict(zip(['b','a'],[seq1(ran2),seq2(ran2)]))
                
                self.hashequal([{1:d1,2:d2}]*3,[{1:d2,2:d1}]*3)
                self.hashnotequal([{1:d1,2:d2}]*3,[{1:d3,2:d1}]*3)
                
    def test_inheritance(self):
        n = 20 # TODO: not scaleable
        superdict = dict(zip([chr(c+ord('a')) for c in range(n)],np.linspace(0,1,n)))
        
        root = parameters.Node(name="root")
        params = root.dtype(random.sample(superdict.items(), n//2))
        root.update(params)
        lst = [params]
        nodes = [root]
        root.reset()
        #print -1
        #print params
        
        n = 2*n//3
        for i in range(n):
            node = nodes[-1].branch(name="Node{}".format(i))
            nodes.append(node)
            
            params = node.dtype(random.sample(superdict.items(), n-i))
            for k,b in zip(params.keys(),np.random.choice([False,True],n-i)):
                if b:
                    params[k] = random.random()

            paramsi = lst[-1].copy()
            paramsi.update(params)
            #print i
            #print params
            #print paramsi
            #print ""
            node.update(params)
            lst.append(paramsi)
        
        for params,node in zip(lst,nodes):
            #print ""
            #print root.tree()
            #print node.name
            #print node.items()
            #print node.todict()
            #print params
            self.assertEqual(params,node)

def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_parameters("test_hash"))
    testSuite.addTest(test_parameters("test_inheritance"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

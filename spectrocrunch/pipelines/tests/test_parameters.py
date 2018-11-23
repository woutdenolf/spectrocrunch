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
from ...utils import hashing
import numpy as np
import random
import collections
import itertools

class Counter(object):

    def __init__(self,start=0):
        self._ctr = start-1
        
    def __call__(self):
        self._ctr+=1
        return self._ctr
        
        
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
    
    def _parameters(self,n):
        letters = [chr(c+ord('a')) for c in range(26)]
        if n>26:
            letters = list(itertools.combinations_with_replacement(letters, (n+26-1)//26))
        letters = letters[:n]
        return dict(zip(letters,np.arange(n)))

    def _parameters_sample(self,superdict,n,overwrite=True):
        params = dict(random.sample(superdict.items(), n))
        if overwrite:
            for k,b in zip(params.keys(),np.random.choice([False,True],n)):
                if b:
                    params[k] = random.random()
        return params
        
    def _generate_parameters(self,n=100):
        superdict = self._parameters(n)
        params_local = [self._parameters_sample(superdict,n//2,overwrite=False)]
        n = 2*n//3
        for i in range(n):
            params_local.append(self._parameters_sample(superdict,n-i))
        return params_local
    
    def _parameters_inherit(self,params_local):
        params_inherit = []
        
        for i,params in enumerate(params_local):
            if params_inherit:
                paramsi = params_inherit[-1].copy()
                paramsi.update(params)
            else:
                paramsi = params.copy()
            params_inherit.append(paramsi)
            
        return params_inherit

    def _generate_nodes_branch(self,params_local):
        nodes = []

        for i,params in enumerate(params_local):
            name = "Node{}".format(i)
            if nodes:
                node = nodes[-1].branch(name=name)
            else:
                node = parameters.Node(name=name)
            nodes.append(node)
            node.update(params)
            
        return nodes

    def _getnodes(self,node):
        nodes = []
        while node is not None:
            nodes.append(node)
            if node.children:
                node = node.children[0]
            else:
                node = None
        return nodes
    
    @staticmethod
    def _lst_move(lst,index_source,before,index_dest):
        if before:
            if index_source>=index_dest:
                pass
            else:
                index_dest -= 1
        else:
            if index_source<=index_dest:
                pass
            else:
                index_dest += 1
        lst.insert(index_dest, lst.pop(index_source))
            
    def _shuffle_nodes(self,nodes):
        # Only works for a single chain of nodes
        params_sorted = [node.todict(full=False) for node in nodes]
        
        for i in range(len(nodes)*2):
            index_source = random.randrange(len(nodes))
            index_dest = random.randrange(len(nodes))
            before = random.choice([True,False])
            # single=True: do not move children (causes forked chain)
            if before:
                _,newnode = nodes[index_dest].insert_before(node=nodes[index_source],single=True)
            else:
                _,newnode = nodes[index_dest].insert_after(node=nodes[index_source],single=True)

            self._lst_move(params_sorted,index_source,before,index_dest)
            self._lst_move(nodes,index_source,before,index_dest)
            
            if newnode is not None:
                nodes.insert(0,newnode)
                params_sorted.insert(0,{})
        
        parameters_inherit = self._parameters_inherit(params_sorted)
        return nodes,parameters_inherit
        
    def test_inheritance(self):
        params_local = self._generate_parameters()
        params_inherited = self._parameters_inherit(params_local)
        nodes = self._generate_nodes_branch(params_local)

        for node,params in zip(nodes,params_inherited):
            self.assertEqual(params,node.todict())
    
    def _assert_nodes_order(self,nodes):
        lst1 = [node.name for node in nodes[0].root.iter_down()]
        lst2 = [node.name for node in nodes]
        self.assertEqual(lst1,lst2)
    
    def test_insert_parameters(self):
        for nlevels in [1,2,4,10]:
            params_local = self._generate_parameters(n=nlevels)
            nodes = self._generate_nodes_branch(params_local)
            nodes,params_inherited = self._shuffle_nodes(nodes)
            self._assert_nodes_order(nodes)
            for node,params in zip(nodes,params_inherited):
                self.assertEqual(params,node)
                
    def _checknodes(self,nodes):
        root = nodes[0].root
        
        lst = list(root.iter_down())
        self.assertEqual(len(lst),len(nodes))
        
        lst = set(id(child) for node in nodes for child in node.children)
        lst.add(id(root))
        self.assertEqual(len(lst),len(nodes))
        
        lst = set(id(node.parent) for node in nodes if id(node)!=id(root))
        lst |= set(id(node) for node in nodes if id(node)!=id(root) and not node.children)
        self.assertEqual(len(lst),len(nodes))
        
    def test_insert_connections(self):
        for nlevels in [1,2,4,10]:
            root = parameters.Node(name="root")
            nodes = [root]
            nodeslast = [root]
            ctr = Counter()
            for i in range(nlevels):
                nodeslast = [node.branch(name="Node{}".format(ctr())) for node in nodeslast for i in range(2)]
                nodes.extend(nodeslast)
            self._checknodes(nodes)

            for j in range(50):
                func = random.choice(["insert_after","insert_before"])
                single = random.choice([True,False])
                node1 = random.choice(nodes)
                node2 = random.choice(nodes)
                _,newnode = getattr(node1,func)(node=node2,single=single)
                if newnode is not None:
                    if not any(newnode is node for node in nodes):
                        nodes.append(newnode)
                self._checknodes(nodes)

    def _random_tree(self,name="Node",nlevels=2):
        root = parameters.Node(name="{}0".format(name))
        nodeslast = [root]
        nodes = [root]
        ctr = Counter(1)
        for i in range(nlevels):
            nodeslast = [node.branch(name="{}{}".format(name,ctr())) for node in nodeslast for i in range(2)]
            nodes.extend(nodeslast)
        
        n = len(nodes)
        superdict = self._parameters(n)
        params = self._parameters_sample(superdict,n//2,overwrite=False)
        nodes[0].update(params)
        n = 2*n//3
        for i,node in enumerate(nodes[1:]):
            params = self._parameters_sample(superdict,max(n-i,1))
            node.update(params)
        
        return root
        
    def debug1(self):
        rootA = self._random_tree(name="A",nlevels=2)
        rootB = self._random_tree(name="B",nlevels=1)

        nodesA = list(rootA.iter_down())
        nodesB = list(rootB.iter_down())
        nodes = []
        for nodeB in nodesB:
            nodeA = random.choice(nodesA)
            nodeA.update(param=nodeB)
            nodes.append(nodeA)
            
        rootA.reset_state(recursive=True)
        rootB.reset_state(recursive=True)
        
        random.choice(nodes)["param"]["change"] = True
        print ""
        print rootA.tree()
        print ""
        print rootB.tree()
    
    def debug2(self):
        root = parameters.Node(name="root")
        node1 = root.branch(name="node1")
        node2 = node1.branch(name="node2")
        
        node = node1
        node["var"] = 10
        node["fit"] = {"linear":True,"Lines":["Fe-K","Ca-K"]}
        #node["fit"]["test"] = {'a':1}
        #root.reset_state(recursive=True)
        #node["fit"]["linear"] = False
        
        node = node2
        node["fit"] = {}
        node["fit"]["linear"] = False
        #node2["fit"]["test"]['a'] = 2
        print root.tree()
        
        #print node1

def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_parameters("test_hash"))
    testSuite.addTest(test_parameters("test_inheritance"))
    testSuite.addTest(test_parameters("test_insert_connections"))
    testSuite.addTest(test_parameters("test_insert_parameters"))
    #testSuite.addTest(test_parameters("debug2"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

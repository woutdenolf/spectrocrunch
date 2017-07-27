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

import unittest

from . import genindexing

from .. import indexing

import numpy as np

import itertools

class test_indexing(unittest.TestCase):

    def test_list(self):
        for k in range(50):
            n = 100
            lst = range(n)
            for index in genindexing.genindexing(n,advanced=False):
                slst = lst[index]
                self.assertEqual(indexing.nonchanging(index,n),slst==lst)
    
    def test_numpy(self):

        for k in range(50):
            n1 = 4
            n2 = 10
            arr = np.zeros((n1,n2))

            index = genindexing.genindexingn((n1,n2),advanced=True)
            for i in index:
                if i[0] is None:
                    i = (None,slice(None),i[0])

                sarr = arr[i]
                self.assertEqual(indexing.nonchanging(i,[n1,n2]),np.array_equal(sarr,arr))

            for index1 in genindexing.genindexing(n1,advanced=True):
                index = index1
                sarr = arr[index]

                self.assertEqual(indexing.nonchanging(index,n1),np.array_equal(sarr,arr))

                for index2 in genindexing.genindexing(n2,advanced=True):
                    index = (index1,index2)
                    if not genindexing.valid(index):
                        continue
                    if index1 is None:
                        index = (None,slice(None),index2)

                    sarr = arr[index]
                    self.assertEqual(indexing.nonchanging(index,[n1,n2]),np.array_equal(sarr,arr))

    def test_replacefull(self):
        a = np.arange(20*30*40).reshape(20,30,40)

        ind = indexing.replacefull(Ellipsis,3,0)
        np.testing.assert_array_equal(a,a[ind])
        ind = indexing.replacefull(Ellipsis,3,1)
        np.testing.assert_array_equal(a,a[ind])
        ind = indexing.replacefull(Ellipsis,3,2)
        np.testing.assert_array_equal(a,a[ind])

        ind = indexing.replacefull((slice(None),Ellipsis),3,-1)
        np.testing.assert_array_equal(a,a[ind])
        ind = indexing.replacefull((Ellipsis,slice(None)),3,2)
        np.testing.assert_array_equal(a,a[ind])
        ind = indexing.replacefull((Ellipsis,slice(None)),3,0)
        np.testing.assert_array_equal(a,a[ind])

        ind = indexing.replacefull((slice(5,7),Ellipsis),3,0)
        np.testing.assert_array_equal(a,a[ind])
        ind = indexing.replacefull((Ellipsis,slice(5,7)),3,2)
        np.testing.assert_array_equal(a,a[ind])
        
        ind = indexing.replacefull((slice(5,7),slice(None),slice(None)),3,0)
        np.testing.assert_array_equal(a,a[ind])
        ind = indexing.replacefull((slice(None),slice(5,7),slice(None)),3,1)
        np.testing.assert_array_equal(a,a[ind])
        ind = indexing.replacefull((slice(None),slice(None),slice(5,7)),3,2)
        np.testing.assert_array_equal(a,a[ind])

    def test_expand(self):
        self.assertEqual(indexing.expand(  (0,slice(0,10))  ,3), (0,slice(0,10),slice(None)) )
        self.assertEqual(indexing.expand(  0  ,3), (0,slice(None),slice(None)) )
        self.assertEqual(indexing.expand(  slice(0,10)  ,3), (slice(0,10),slice(None),slice(None)) )
        self.assertEqual(indexing.expand(  (None,slice(0,10))  ,3), (None,slice(0,10),slice(None),slice(None)) )
        self.assertEqual(indexing.expand(  (slice(0,10),None)  ,3), (slice(0,10),None,slice(None),slice(None)) )
        self.assertEqual(indexing.expand(  None  ,3), (None,slice(None),slice(None),slice(None)) )
        self.assertEqual(indexing.expand(  (None,None)  ,3), (None,None,slice(None),slice(None),slice(None)) )
        self.assertEqual(indexing.expand(  (None,Ellipsis,None)  ,3), (None,slice(None),slice(None),slice(None),None) )
        self.assertEqual(indexing.expand(  Ellipsis  ,3), (slice(None),slice(None),slice(None)) )

    def test_axisindex(self):
        ndim = 10

        ind = (0,Ellipsis,slice(0,5,5),None,2,slice(1,2),None)
        ind = indexing.expand(ind,ndim)

        self.assertEqual(ind[indexing.axisindex(ind,-1,ndim)],slice(1,2))
        self.assertEqual(ind[indexing.axisindex(ind,8,ndim)],2)
        self.assertEqual(ind[indexing.axisindex(ind,7,ndim)],slice(0,5,5))
        self.assertEqual(ind[indexing.axisindex(ind,1,ndim)],slice(None))

    def test_extract_dimnonchanging(self):
        a = np.arange(10*20*30*40).reshape((10,20,30,40))

        index = (slice(1,2),range(5),slice(3,4),range(5))
        index1 = indexing.extract_dimchanging(index)
        index2,iadv = indexing.extract_dimnonchanging(index)
        self.assertEqual(iadv,0)
        np.testing.assert_array_equal(a[index],a[index1][index2])

        index = (slice(1,2),range(5),range(5),slice(3,4))
        index1 = indexing.extract_dimchanging(index)
        index2,iadv = indexing.extract_dimnonchanging(index)
        self.assertEqual(iadv,1)
        np.testing.assert_array_equal(a[index],a[index1][index2])

        index = (slice(1,2),range(5),0,range(5))
        index1 = indexing.extract_dimchanging(index)
        index2,iadv = indexing.extract_dimnonchanging(index)
        self.assertEqual(iadv,1)
        np.testing.assert_array_equal(a[index],a[index1][index2])

        index = (range(5),0,range(5),slice(1,2))
        index1 = indexing.extract_dimchanging(index)
        index2,iadv = indexing.extract_dimnonchanging(index)
        self.assertEqual(iadv,0)
        np.testing.assert_array_equal(a[index],a[index1][index2])

        index = (0,0,slice(1,2),range(5))
        index1 = indexing.extract_dimchanging(index)
        index2,iadv = indexing.extract_dimnonchanging(index)
        self.assertEqual(iadv,0)
        np.testing.assert_array_equal(a[index],a[index1][index2])

        index = (slice(None),np.newaxis,range(5),range(5),0,np.newaxis)
        index1 = indexing.extract_dimchanging(index)
        index2,iadv = indexing.extract_dimnonchanging(index)
        self.assertEqual(iadv,2)
        np.testing.assert_array_equal(a[index],a[index1][index2])

        index = (np.newaxis,range(5),0,range(5),0,np.newaxis)
        index1 = indexing.extract_dimchanging(index)
        index2,iadv = indexing.extract_dimnonchanging(index)
        self.assertEqual(iadv,1)
        np.testing.assert_array_equal(a[index],a[index1][index2])

        index = (slice(None),np.newaxis,0,0,np.newaxis)
        index1 = indexing.extract_dimchanging(index)
        index2,iadv = indexing.extract_dimnonchanging(index)
        self.assertEqual(iadv,None)
        np.testing.assert_array_equal(a[index],a[index1][index2])

        ndim = 4
        for nsingleton in range(ndim):
            for nlist in range(ndim-nsingleton):
                for nother in range(ndim-nsingleton-nlist):
                    for nnew in range(0,2):
                        p = [range(5)]*nlist + [0]*nsingleton + [slice(None)]*nother + [np.newaxis]*nnew
                        for index in itertools.permutations(p):
                            #index = indexing.expand(index,ndim=ndim)
                            index1 = indexing.extract_dimchanging(index)
                            index2,iadv = indexing.extract_dimnonchanging(index)
                            np.testing.assert_array_equal(a[index],a[index1][index2])

    def test_extract_newaxis(self):
        a = np.arange(10*20*30*40).reshape((10,20,30,40))

        ndim = 4
        for nsingleton in range(ndim):
            for nlist in range(ndim-nsingleton):
                for nother in range(ndim-nsingleton-nlist):
                    for nnew in range(0,2):
                        p = [range(5)]*nlist + [0]*nsingleton + [slice(None)]*nother + [np.newaxis]*nnew
                        for index1 in itertools.permutations(p):
                            index1 = indexing.expand(index1,ndim=ndim)
                            index2,o = indexing.extract_newaxis(index1)
                            np.testing.assert_array_equal(a[index1],o(a[index2]))


    def _slicegen(self,args):
        return np.arange(np.prod(args[0])).reshape(args[0])+args[1]

    def test_funcslicing(self):
        args = [[(20,30,40),i*1000] for i in range(10)]

        ndim = 4
        for nsingleton in range(ndim):
            for nlist in range(ndim-nsingleton):
                for nother in range(ndim-nsingleton-nlist):
                    for nnew in range(0,2):
                        p = [range(5)]*nlist + [0]*nsingleton + [slice(None)]*nother + [np.newaxis]*nnew
                        for index in itertools.permutations(p):
                            for axis in range(-ndim,ndim):
                                data = np.stack(map(self._slicegen,args),axis=axis)
                                shapefull = data.shape

                                data2 = indexing.slicedstack(self._slicegen,args,index,ndim,shapefull=shapefull,axis=axis)
                                np.testing.assert_array_equal(data[index],data2)

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_indexing("test_list"))
    testSuite.addTest(test_indexing("test_numpy"))
    testSuite.addTest(test_indexing("test_replacefull"))
    testSuite.addTest(test_indexing("test_expand"))
    testSuite.addTest(test_indexing("test_axisindex"))
    testSuite.addTest(test_indexing("test_extract_dimnonchanging"))
    testSuite.addTest(test_indexing("test_extract_newaxis"))
    testSuite.addTest(test_indexing("test_funcslicing"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

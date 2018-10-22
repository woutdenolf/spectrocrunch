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
from .. import listtools

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
                if i[0] is np.newaxis:
                    i = (np.newaxis,slice(None),i[0])
                sarr = arr[i]
                self.assertEqual(indexing.nonchanging(i,[n1,n2]),np.array_equal(sarr,arr))

            for index1 in genindexing.genindexing(n1,advanced=True):
                index = index1
                sarr = arr[index]

                self.assertEqual(indexing.nonchanging(index,n1),np.array_equal(sarr,arr))

                for index2 in genindexing.genindexing(n2,advanced=True):
                    index = (index1,index2)
                    if index1 is np.newaxis:
                        index = (np.newaxis,slice(None),index2)
                    if not genindexing.valid(index,(n1,n2)):
                        continue
                    sarr = arr[index]
                    self.assertEqual(indexing.nonchanging(index,[n1,n2]),np.array_equal(sarr,arr))

    def test_replacefull(self):
        a = np.arange(20*30*40).reshape(20,30,40)

        ind = indexing.replacefull(Ellipsis,3,[0])
        np.testing.assert_array_equal(a,a[ind])
        ind = indexing.replacefull(Ellipsis,3,[1])
        np.testing.assert_array_equal(a,a[ind])
        ind = indexing.replacefull(Ellipsis,3,[2])
        np.testing.assert_array_equal(a,a[ind])

        ind = indexing.replacefull((slice(None),Ellipsis),3,[-1])
        np.testing.assert_array_equal(a,a[ind])
        ind = indexing.replacefull((Ellipsis,slice(None)),3,[2])
        np.testing.assert_array_equal(a,a[ind])
        ind = indexing.replacefull((Ellipsis,slice(None)),3,[0])
        np.testing.assert_array_equal(a,a[ind])

        ind = indexing.replacefull((slice(5,7),Ellipsis),3,[0])
        np.testing.assert_array_equal(a,a[ind])
        ind = indexing.replacefull((Ellipsis,slice(5,7)),3,[2])
        np.testing.assert_array_equal(a,a[ind])
        
        ind = indexing.replacefull((slice(5,7),slice(None),slice(None)),3,[0])
        np.testing.assert_array_equal(a,a[ind])
        ind = indexing.replacefull((slice(None),slice(5,7),slice(None)),3,[1])
        np.testing.assert_array_equal(a,a[ind])
        ind = indexing.replacefull((slice(None),slice(None),slice(5,7)),3,[2])
        np.testing.assert_array_equal(a,a[ind])

    def test_expand(self):
        self.assertEqual(indexing.expand(  (0,slice(0,10))  ,3), (0,slice(0,10),slice(None)) )
        self.assertEqual(indexing.expand(  0  ,3), (0,slice(None),slice(None)) )
        self.assertEqual(indexing.expand(  slice(0,10)  ,3), (slice(0,10),slice(None),slice(None)) )
        self.assertEqual(indexing.expand(  (np.newaxis,slice(0,10))  ,3), (np.newaxis,slice(0,10),slice(None),slice(None)) )
        self.assertEqual(indexing.expand(  (slice(0,10),np.newaxis)  ,3), (slice(0,10),None,slice(None),slice(None)) )
        self.assertEqual(indexing.expand(  np.newaxis  ,3), (np.newaxis,slice(None),slice(None),slice(None)) )
        self.assertEqual(indexing.expand(  (np.newaxis,np.newaxis)  ,3), (np.newaxis,np.newaxis,slice(None),slice(None),slice(None)) )
        self.assertEqual(indexing.expand(  (np.newaxis,Ellipsis,np.newaxis)  ,3), (np.newaxis,slice(None),slice(None),slice(None),np.newaxis) )
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
                        for r in [0,5]:
                            p = [range(r)]*nlist + [0]*nsingleton + [slice(None)]*nother + [np.newaxis]*nnew
                            for index in itertools.permutations(p):
                                #index = indexing.expand(index,ndim=ndim) # apparently not needed?
                                index1 = indexing.extract_dimchanging(index)
                                index2,iadv = indexing.extract_dimnonchanging(index)
                                np.testing.assert_array_equal(a[index],a[index1][index2])

    def test_nonchangingdims(self):
        ndim = 4
        axes = [0,1]
        shape = (4,5,6,7)

        self.assertTrue(indexing.nonchangingdims((slice(None),slice(None),slice(None),slice(None)),ndim,axes,shape=shape))
        self.assertTrue(indexing.nonchangingdims((slice(None),slice(None),[0,0,0],[0,0,0]),ndim,axes,shape=shape))
        self.assertTrue(indexing.nonchangingdims((slice(None),slice(None),0,[0,0,0]),ndim,axes,shape=shape))
        self.assertFalse(indexing.nonchangingdims((slice(None),slice(None),0,np.newaxis,[0,0,0]),ndim,axes,shape=shape))
        self.assertTrue(indexing.nonchangingdims(([True]*4,slice(None),slice(None),slice(None)),ndim,axes,shape=shape))
        self.assertTrue(indexing.nonchangingdims((range(4),slice(None),slice(None),slice(None)),ndim,axes,shape=shape))
        self.assertFalse(indexing.nonchangingdims(([False,True,True,True],slice(None),slice(None),slice(None)),ndim,axes,shape=shape))
        self.assertFalse(indexing.nonchangingdims(([1,0,2,3],slice(None),slice(None),slice(None)),ndim,axes,shape=shape))

    def test_extract_newaxis(self):
        a = np.arange(10*20*30*40).reshape((10,20,30,40))
        self.execute_test_extract_newaxis(a,[0,5])

        a = np.ones((1,1,1,1))
        self.execute_test_extract_newaxis(a,[0])

    def execute_test_extract_newaxis(self,a,rrange):
        ndim = a.ndim
        for nsingleton in range(ndim):
            for nlist in range(ndim-nsingleton):
                for nother in range(ndim-nsingleton-nlist):
                    for nnew in range(0,2):
                        for r in rrange:
                            p = [range(r)]*nlist + [0]*nsingleton + [slice(None)]*nother + [np.newaxis]*nnew
                            for index1 in itertools.permutations(p):
                                index1 = indexing.expand(index1,ndim=ndim)
                                index2,o = indexing.extract_newaxis(index1)
                                np.testing.assert_array_equal(a[index1],o(a[index2]))
                         
    def test_shape_afterindexing(self):
        a = np.arange(10*20*30*40).reshape((10,20,30,40))
        self.try_shape_afterindexing(a,[0,5])

        a = np.ones((1,1,1,1))
        self.try_shape_afterindexing(a,[0])

    def try_shape_afterindexing(self,a,rrange):
        s1 = list(a.shape) # original shape
        ndim = a.ndim
        for nsingleton in range(ndim):
            for nlist in range(ndim-nsingleton):
                for nother in range(ndim-nsingleton-nlist):
                    for nnew in range(0,2):
                        for r in rrange:
                            p = [range(r)]*nlist + [0]*nsingleton + [slice(None)]*nother + [np.newaxis]*nnew
                            for index in itertools.permutations(p):
                                # Shape after indexing
                                s3 = a[index].shape

                                # Expected shape after indexing
                                s2 = indexing.shape_afterindexing(s1,index)
                                self.assertEqual(s2,s3)

    def test_replacefull_transform(self):
        a = np.arange(10*20*30*40).reshape((10,20,30,40))
        ndim = a.ndim

        for nsingleton in range(ndim):
            for nlist in range(ndim-nsingleton):
                for nother in range(ndim-nsingleton-nlist):
                    for nnew in range(0,1):
                        for r in [5]:
                            p = [range(r)]*nlist + [0]*nsingleton + [slice(None)]*nother + [np.newaxis]*nnew
                            for index in itertools.permutations(p):
                                for n in range(1,ndim):
                                    
                                    for fullaxes in itertools.combinations(range(ndim),n):
                                        indexfull,postindexfull,singletonindex = indexing.replacefull_transform(index,fullaxes,ndim,restoreadvanced=False)
 
                                        ind1 = indexing.expand(index,ndim)
                                        ind2 = indexfull
                                        ind3 = [i for i in ind1 if i is not np.newaxis]
                                        b1 = [isinstance(i, list) for i in ind1]
                                        b2 = [not isinstance(i, list) for i in ind2]
                                        bseladv = np.logical_and(b1,b2)

                                        if any(bseladv):
                                            selind = [[ i if isinstance(ind3[axis],list) else 0 for axis in fullaxes] for i in range(r)]
                                        else:
                                            selind = [[0]*n]
         
                                        # Normal indexing
                                        b = a[index]
                                        
                                        # Indexing without fullaxes (fullaxes might be increased because of advanced indexing)
                                        c = postindexfull(a[indexfull])
                                        
                                        # Index fullaxes after indexing without fullaxes
                                        d = singletonindex(c,selind)
                                        
                                        if any(bseladv):
                                            if r in d[0].shape:
                                                i = d[0].shape.index(r)
                                                singletonindex = indexing.op_singletonindex([i],[False])
                                                d = [singletonindex(x,[[i]])[0] for i,x in enumerate(d)]
                                            d = np.stack(d,axis=b.shape.index(r))
                                        else:
                                            d = d[0]
                                        
                                        self.assertTrue(b.ndim,d.ndim)

                                        # Index fullaxes after normal indexing
                                        selaxes = np.logical_and( np.asarray(d.shape)==1,np.asarray(b.shape)!=1 )
                                        if any(selaxes):
                                            selaxes = np.arange(b.ndim)[selaxes]
                                            restore = [True]*len(selaxes)
                                            ind = [0]*len(selaxes)
                                                
                                            singletonindex = indexing.op_singletonindex(selaxes,restore)
                                            b = singletonindex(b,[ind])[0]

                                        np.testing.assert_array_equal(b,d)

    def _slicegen(self,args):
        return np.arange(np.prod(args[0])).reshape(args[0])+args[1]

    def test_getitem(self):
        args = [[(10,11,12),i*1000] for i in range(8)]
        self.execute_test_getitem(self._slicegen,args,[0,5])
        
        args = [[(1,1,1),i*1000] for i in range(1)]
        self.execute_test_getitem(self._slicegen,args,[0])

        args = [[(10,11,12),i*1000] for i in range(8)]
        special = [(([0], slice(None, None, None), slice(None, None, None), [-1]),-1),\
                   ((0,0,0,[0]),-1),\
                   ((0, slice(None, None, None), slice(None, None, None), [0]),0)]

        ndim = 4
        for index,axis in special:
            data = np.stack(map(self._slicegen,args),axis=axis)
            shapefull = data.shape
            data2 = indexing.getitem(self._slicegen,args,index,ndim,shapefull=shapefull,axis=axis)
            np.testing.assert_array_equal(data[index],data2)

    def test_setitem(self):
        args = [[(10,11,12),i*1000] for i in range(8)]
        self.execute_test_setitem(self._slicegen,args,[1,5])
        
        args = [[(10,11,12),i*1000] for i in range(8)]
        special = [((slice(0, 2, None), slice(2, 3, None), -1, slice(None, None, None)),-1),
                   ((0,0,0,0),-1),
                   ((0,[0,1],0,0),-1)]

        ndim = 4
        for index,axis in special:
            data = np.stack(map(self._slicegen,args),axis=axis)
            shapefull = data.shape
            if axis<0:
                j = axis+ndim
            else:
                j = axis
            selector = lambda i: data[(slice(None),)*j+(i,)+(slice(None),)*(ndim-j-1)]
            args2 = range(shapefull[axis])
            keep = data[index]
            new = np.max(keep) + np.arange(keep.size).reshape(keep.shape)
            indexing.setitem(selector,args2,index,ndim,new,shapefull=shapefull,axis=axis)
            np.testing.assert_array_equal(data[index],new)
            
            
    def execute_test_getitem(self,generator,args,rrange):
        ndim = 4
        for nsingleton in range(ndim):
            for nlist in range(ndim-nsingleton):
                for nother in range(ndim-nsingleton-nlist):
                    for nnew in range(0,2):
                        for r in rrange:
                            p = [range(r)]*nlist + [0]*nsingleton + [slice(None)]*nother + [np.newaxis]*nnew
                            for index in itertools.permutations(p):
                                for axis in range(-ndim,ndim):
                                    data = np.stack(map(generator,args),axis=axis)
                                    shapefull = data.shape
                                    data2 = indexing.getitem(generator,args,index,ndim,shapefull=shapefull,axis=axis)
                                    np.testing.assert_array_equal(data[index],data2)

    def execute_test_setitem(self,generator,args,rrange):
        ndim = 4
        for nsingleton in range(ndim):
            for nlist in range(ndim-nsingleton):
                for nother in range(ndim-nsingleton-nlist):
                    for nnew in [0]:
                        for r in rrange:
                            p = [range(r)]*nlist + [0]*nsingleton + [slice(None)]*nother + [np.newaxis]*nnew
                            for index in itertools.permutations(p):
                                for axis in range(-ndim,ndim):
                                    data = np.stack(map(generator,args),axis=axis)
                                    shapefull = data.shape
                                    if axis<0:
                                        j = axis+ndim
                                    else:
                                        j = axis
                                    selector = lambda i: data[(slice(None),)*j+(i,)+(slice(None),)*(ndim-j-1)]
                                    args2 = range(shapefull[axis])
                                    keep = data[index]
                                    new = np.max(keep) + np.arange(keep.size).reshape(keep.shape)
                                    indexing.setitem(selector,args2,index,ndim,new,shapefull=shapefull,axis=axis)
                                    np.testing.assert_array_equal(data[index],new)
                                    
                                    
def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_indexing("test_list"))
    testSuite.addTest(test_indexing("test_numpy"))
    testSuite.addTest(test_indexing("test_replacefull"))
    testSuite.addTest(test_indexing("test_expand"))
    testSuite.addTest(test_indexing("test_axisindex"))
    testSuite.addTest(test_indexing("test_nonchangingdims"))
    testSuite.addTest(test_indexing("test_extract_dimnonchanging"))
    testSuite.addTest(test_indexing("test_extract_newaxis"))
    testSuite.addTest(test_indexing("test_shape_afterindexing"))
    testSuite.addTest(test_indexing("test_replacefull_transform"))
    testSuite.addTest(test_indexing("test_getitem"))
    testSuite.addTest(test_indexing("test_setitem"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

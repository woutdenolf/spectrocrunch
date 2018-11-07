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
import numpy as np

from .. import regulargrid
from .. import axis
from ...utils import units
from ...utils.tests import genindexing

class test_regulargrid(unittest.TestCase):

    def _generate(self,stackaxis,shape=(2,3,4),stackdim=0):
        dtype = np.float32
        stackaxis = axis.factory(stackaxis,type='ordinal')
        axes = [axis.factory(range(n)) for n in shape]
        axes.insert(stackdim,stackaxis)
        shape = [len(ax) for ax in axes]
        data = np.random.normal(size=shape).astype(dtype)*100
        return regulargrid.RegularGrid(axes,data,stackdim=stackdim)
    
    def test_indexing(self):
        stackaxis = ['Fe-K','Si-K','Al-K','S-K','Ce-L']
        for i in range(3):
            grid = self._generate(stackaxis,shape=(6,8),stackdim=i)
            self._check_grid(grid)
        
    def test_interpolate(self):
        for degree in [0,1,2]:
            rtol = 1e-3
            
            # 1D
            stackaxis = ['Fe-K','Si-K','Al-K','S-K','Ce-L']
            grid = self._generate(stackaxis,shape=(),stackdim=0)
            data = grid.interpolate('Fe-K')
            np.testing.assert_allclose(grid.data[0],data,rtol=rtol)
            data = grid.interpolate(['Fe-K','Ce-L'],degree=degree)
            np.testing.assert_allclose(grid.data[[0,-1]],data,rtol=rtol)
            
            # 2D
            grid = self._generate(stackaxis,shape=(6,),stackdim=0)
            data = grid.interpolate('Fe-K',None,degree=degree)
            np.testing.assert_allclose(grid.data[np.newaxis,0,:],data,rtol=rtol)
            for asgrid in [True,False]:
                data = grid.interpolate(['Fe-K','S-K','Ce-L'],[0,2,3],degree=degree,asgrid=asgrid)
                ind = [0,3,4],[0,2,3]
                if asgrid:
                    ind = np.meshgrid(*ind, indexing='ij')
                np.testing.assert_allclose(grid.data[ind],data,rtol=rtol)
            
            # 3D
            grid = self._generate(stackaxis,shape=(6,8),stackdim=1)
            data = grid.interpolate(None,'Fe-K',None,degree=degree)
            np.testing.assert_allclose(grid.data[:,np.newaxis,0,:],data,rtol=rtol)
            for asgrid in [True,False]:
                data = grid.interpolate([1,3,4],['Fe-K','S-K','Ce-L'],[0,2,3],degree=degree,asgrid=asgrid)
                ind = [1,3,4],[0,3,4],[0,2,3]
                if asgrid:
                    ind = np.meshgrid(*ind, indexing='ij')
                np.testing.assert_allclose(grid.data[ind],data,rtol=rtol)
                
    def _check_grid(self,grid):
        data = grid.values
        self.assertEqual(grid.shape,data.shape)
        self.assertEqual(grid.ndim,data.ndim)
        self.assertEqual(grid.size,data.size)
        np.testing.assert_array_equal(grid[:],data)

        indices = genindexing.genindexingn(data.shape,advanced=False,eco=False,nmax=50)
        for index in indices:
            np.testing.assert_array_equal(grid[index],data[index])
 
        for a,b in zip(grid,data):
            np.testing.assert_array_equal(a,b)
            
def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_regulargrid("test_indexing"))
    testSuite.addTest(test_regulargrid("test_interpolate"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

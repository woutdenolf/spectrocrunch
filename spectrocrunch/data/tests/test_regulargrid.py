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
from testfixtures import TempDirectory
import os
import contextlib

from .. import regulargrid
from ...utils import units
from ...io import nxfs
from ...utils.tests import genindexing
from .. import nxprocess

class test_regulargrid(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()
    
    @contextlib.contextmanager
    def _nxprocess(self,method=None):
        h5filename = os.path.join(self.dir.path,'test.h5')
        root = nxfs.Path('/',h5file=h5filename).nxroot()
        entry = root.new_nxentry()
        process,new = entry.nxprocess('fromraw',parameters={'a':1,'b':2},previous=None)
        info = {}
        
        shape = (2,10,13)
        self.stackdim = 0
        positioners = entry.positioners()
        positioners.add_axis('z',range(shape[0]))
        positioners.add_axis('y',range(shape[1]),units='um',title='vertical')
        positioners.add_axis('x',range(shape[2]))
        
        dtype = np.float32
        signals = ['Fe-K','Si-K','Al-K','S-K','Ce-L']
        for detector in range(2):
            detector = process.results.nxdata('detector{:02d}'.format(detector)).mkdir()
            for name in signals:
                data = np.random.normal(size=shape).astype(dtype)
                if method=='crop':
                    data[:,0,:] = np.nan
                    data[:,-2:,:] = np.nan
                    data[:,:,0:2] = np.nan
                    data[:,:,-1] = np.nan
                    info['y_crop'] = positioners.get_axis('y')[1:-2]
                    info['x_crop'] = positioners.get_axis('x')[2:-1]
                detector.add_signal(name,data=data)
            detector.set_axes('z','y','x')

        try:
            yield process,info
        finally:
            #root.remove(recursive=True)
            pass
                
    def test_nxprocess(self):
        with self._nxprocess() as proc:
            proc,info = proc
            grid = regulargrid.NXRegularGrid(proc)
            self._check_grid(grid)
    
    def test_nxdata(self):
        with self._nxprocess() as proc:
            proc,info = proc
            nxdata = proc.results['detector00']
            grid = regulargrid.NXSignalRegularGrid(nxdata.signal)
            self._check_grid(grid)
    
    def test_crop(self):
        with self._nxprocess(method='crop') as proc1:
            proc1,info = proc1
            parameters = {'method':'crop','reference':'Al-K','nanval':np.nan,'sliced':False,'default':'Si-K'}
            proc2 = nxprocess.single_regulargrid(proc1,parameters)
            parameters['sliced'] = True
            parameters['name'] = 'crop2'
            proc3 = nxprocess.single_regulargrid(proc1,parameters)
            
            grid2 = regulargrid.NXRegularGrid(proc2)
            grid3 = regulargrid.NXRegularGrid(proc3)
            np.testing.assert_array_equal(grid2.values,grid3.values)
            for k,v in info.items():
                for ax in grid2.axes:
                    if ax.name==k:
                        np.testing.assert_array_equal(ax.magnitude,v)
                        break
                else:
                    assert(False)

            self.assertEqual(proc1.default.signal.name,parameters['default'])
            
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
    testSuite.addTest(test_regulargrid("test_nxprocess"))
    testSuite.addTest(test_regulargrid("test_nxdata"))
    testSuite.addTest(test_regulargrid("test_crop"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

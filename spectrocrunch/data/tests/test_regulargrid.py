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

class test_regulargrid(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()
    
    def tearDown(self):
        self.dir.cleanup()
    
    @contextlib.contextmanager
    def _nxprocess(self):
        h5filename = os.path.join(self.dir.path,'test.h5')
        root = nxfs.Path('/',h5file=h5filename).nxroot()
        entry = root.new_nxentry()
        process = entry.nxprocess('test',parameters={'a':1,'b':2},previous=None)

        shape = (4,7)
        positioners = entry.positioners()
        positioners.add_axis('y',range(shape[0]),units='um',title='vertical')
        positioners.add_axis('x',range(shape[1]))
        
        dtype = float
        signals = ['Fe-K','Si-K','Al-K','S-K','Ce-L']
        for detector in range(2):
            detector = process.results.nxdata('detector{:02d}'.format(detector)).mkdir()
            for name in signals:
                data = np.random.normal(size=shape).astype(dtype)
                detector.add_signal(name,data=data)
            detector.set_axes('y','x')

        try:
            yield process
        finally:
            root.remove(recursive=True)
            
    def test_nxprocess(self):
        with self._nxprocess() as nxprocess:
            grid = regulargrid.NXRegularGrid(nxprocess)
            self._check_grid(grid)
    
    def test_nxdata(self):
        with self._nxprocess() as nxprocess:
            nxdata = nxprocess.results['detector00']
            grid = regulargrid.NXSignalRegularGrid(nxdata,nxdata.signal)
            self._check_grid(grid)
        
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
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

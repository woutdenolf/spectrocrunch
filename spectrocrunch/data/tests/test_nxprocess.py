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
from ...io import nxfs
from ...utils.tests import genindexing
from .. import nxtask
from .. import axis
from ...utils import units

class test_nxprocess(unittest.TestCase):

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
        yencres = 2.
        xencres = 4.
        yencoff = -4.
        xencoff = -8.

        posmap = np.arange(shape[1])[:,np.newaxis]+np.arange(shape[2])[np.newaxis,:]
        yenc = np.stack([posmap*yencres+yencoff]*shape[self.stackdim],axis=self.stackdim)
        xenc = np.stack([posmap*xencres+xencoff]*shape[self.stackdim],axis=self.stackdim)

        dtype = np.float32
        signals = ['Fe-K','Si-K','Al-K','S-K','Ce-L']
        
        if method=='replace':
            index = tuple([np.random.randint(0,shape[i],10).tolist() for i in range(3)])
            indexinfo = list(index)
            indexinfo.insert(self.stackdim,slice(None))
            info['index'] = tuple(indexinfo)
        elif method=='align':
            info['axes'] = axis.factory(range(shape[0])),\
                           axis.factory(units.Quantity(range(-shape[0]+1,shape[1]),units='um')),\
                           axis.factory(range(-shape[0]+1,shape[2]))
        elif method=='expression':
            info['expression'] = '{}/{arr_iodet}'
            info['skip'] = ['arr_iodet']
            info['select'] = process.results['counters']['arr_iodet']
        elif method=='resample':
            info = {'y':{'counter':'arr_samy','resolution':yencres,'offset':yencoff},
                    'x':{'counter':'arr_samx','resolution':xencres,'offset':xencoff}}
            
        groups = {}
        for group in range(2):
            groups['detector{:02d}'.format(group)] = signals
        groups['counters'] = ['arr_iodet','arr_idet','arr_samy','arr_samx']

        for group,signals in groups.items():
            group = process.results.nxdata(group).mkdir()
            for name in signals:
                if name=='arr_samy':
                    data = yenc
                elif name=='arr_samx':
                    data = xenc
                else:
                    data = np.random.normal(size=shape)
                data = data.astype(dtype)
                if method=='crop':
                    data[:,0,:] = np.nan
                    data[:,-2:,:] = np.nan
                    data[:,:,0:2] = np.nan
                    data[:,:,-1] = np.nan
                    info['y_crop'] = positioners.get_axis('y')[1:-2]
                    info['x_crop'] = positioners.get_axis('x')[2:-1]
                elif method=='replace':
                    data[index] = -1
                elif method=='minlog':
                    data -= np.min(data)*1.1
                elif method=='align':
                    hot = np.max(data)*1.1
                    for i in range(shape[0]):
                        data[i,i,i] = hot
                group.add_signal(name,data=data)
            group.set_axes('z','y','x')

        try:
            yield process,info
        finally:
            #root.remove(recursive=True)
            pass
                
    def test_grid(self):
        with self._nxprocess() as proc:
            proc,info = proc
            grid = regulargrid.NXRegularGrid(proc)
            self._check_grid(grid)
    
            nxdata = proc.results['detector00']
            grid = regulargrid.NXSignalRegularGrid(nxdata.signal)
            self._check_grid(grid)
    
    def test_crop(self):
        with self._nxprocess(method='crop') as proc1:
            proc1,info = proc1
            parameters = {'method':'crop','default':'Si-K','sliced':False,
                          'reference':'Al-K','nanval':np.nan}
            task = nxtask.newtask(parameters,proc1)
            proc2 = task.run()
            proc3 = task.run()
            self.assertEqual(proc2,proc3)
            
            parameters['sliced'] = True
            parameters['name'] = 'crop2'
            task = nxtask.newtask(parameters,proc1)
            proc4 = task.run()
            self.assertNotEqual(proc2,proc4)
            
            grid1 = regulargrid.NXRegularGrid(proc1)
            grid2 = regulargrid.NXRegularGrid(proc2)
            grid4 = regulargrid.NXRegularGrid(proc4)
            self.assertEqual(set([sig.name for sig in grid1.signals]),
                             set([sig.name for sig in grid2.signals]))
            self.assertFalse(np.isnan(grid2.values).any())
            np.testing.assert_array_equal(grid2.values,grid4.values)
            for k,v in info.items():
                for ax in grid2.axes:
                    if ax.name==k:
                        np.testing.assert_array_equal(ax.magnitude,v)
                        break
                else:
                    assert(False)

            self.assertEqual(proc1.default.signal.name,parameters['default'])
    
    def test_replace(self):
        with self._nxprocess(method='replace') as proc1:
            proc1,info = proc1
            parameters = {'method':'replace','default':'Si-K','sliced':False,
                          'old':-1,'new':-2}
            task = nxtask.newtask(parameters,proc1)
            proc2 = task.run()
            proc3 = task.run()
            self.assertEqual(proc2,proc3)
            
            parameters['sliced'] = True
            parameters['name'] = 'replace2'
            task = nxtask.newtask(parameters,proc1)
            proc4 = task.run()
            self.assertNotEqual(proc2,proc4)
            
            grid1 = regulargrid.NXRegularGrid(proc1)
            grid2 = regulargrid.NXRegularGrid(proc2)
            grid4 = regulargrid.NXRegularGrid(proc4)
            self.assertEqual(set([sig.name for sig in grid1.signals]),
                             set([sig.name for sig in grid2.signals]))
            np.testing.assert_array_equal(grid2.values,grid4.values)
            np.testing.assert_array_equal(grid1.values[info['index']],parameters['old'])
            np.testing.assert_array_equal(grid2.values[info['index']],parameters['new'])

            self.assertEqual(proc1.default.signal.name,parameters['default'])
    
    def test_minlog(self):
        with self._nxprocess(method='minlog') as proc1:
            proc1,info = proc1
            parameters = {'method':'minlog','sliced':False}
            task = nxtask.newtask(parameters,proc1)
            proc2 = task.run()
            proc3 = task.run()
            self.assertEqual(proc2,proc3)
            
            parameters['sliced'] = True
            parameters['name'] = 'minlog2'
            task = nxtask.newtask(parameters,proc1)
            proc4 = task.run()
            self.assertNotEqual(proc2,proc4)
            
            grid1 = regulargrid.NXRegularGrid(proc1)
            grid2 = regulargrid.NXRegularGrid(proc2)
            grid4 = regulargrid.NXRegularGrid(proc4)
            
            self.assertEqual(set([sig.name for sig in grid1.signals]),
                             set([sig.name for sig in grid2.signals]))
            np.testing.assert_array_equal(grid2.values,grid4.values)
            np.testing.assert_array_equal(-np.log(grid1),grid4.values)
            
    def test_align(self):
        with self._nxprocess(method='align') as proc1:
            proc1,info = proc1
            parameters = {'method':'align','alignmethod':'max','reference':'Fe-K','refimageindex':0}
            task = nxtask.newtask(parameters,proc1)
            proc2 = task.run()
            proc3 = task.run()
            self.assertEqual(proc2,proc3)
            
            grid2 = regulargrid.NXRegularGrid(proc2)
            axes = grid2.axes
            axes.pop(grid2.stackdim)
            for ax1,ax2 in zip(info['axes'],axes):
                self.assertEqual(ax1,ax2)
    
    def test_expression(self):
        with self._nxprocess(method='expression') as proc1:
            proc1,info = proc1
            parameters = {'method':'expression','expression':info['expression'],'skip':info['skip'],'sliced':False}
            task = nxtask.newtask(parameters,proc1)
            proc2 = task.run()
            proc3 = task.run()
            self.assertEqual(proc2,proc3)
            
            parameters['sliced'] = True
            parameters['name'] = 'expression2'
            task = nxtask.newtask(parameters,proc1)
            proc4 = task.run()
            self.assertNotEqual(proc2,proc4)

            grid1 = regulargrid.NXRegularGrid(proc1)
            grid2 = regulargrid.NXRegularGrid(proc2)
            self._check_axes(grid1,grid2)

            index = grid1.locate(info['select'],None,None,None)
            norm = grid1[index]
            inorm = index[grid1.stackdim]
            for i in range(grid1.shape[grid1.stackdim]):
                index = list(index)
                index[grid1.stackdim] = i
                index = tuple(index)
                data = grid1[index]
                if grid1.signals[i].name not in info['skip']:
                    data = data/norm
                np.testing.assert_array_equal(data,grid2[index])
    
    def test_resample(self):
        with self._nxprocess(method='resample') as proc1:
            proc1,info = proc1
            parameters = {'method':'resample','encoders':info}
            task = nxtask.newtask(parameters,proc1)
            proc2 = task.run()
            proc3 = task.run()
            self.assertEqual(proc2,proc3)

    def _check_axes(self,grid1,grid2):
        for ax1,ax2 in zip(grid1.axes,grid2.axes):
            if ax1.type=='quantitative':
                self.assertEqual(ax1,ax2)
            else:
                self.assertEqual(len(ax1),len(ax2))
                for v1,v2 in zip(ax1,ax2):
                    self.assertEqual(v1.name,v2.name)
                    
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
    testSuite.addTest(test_nxprocess("test_grid"))
    testSuite.addTest(test_nxprocess("test_crop"))
    testSuite.addTest(test_nxprocess("test_replace"))
    testSuite.addTest(test_nxprocess("test_minlog"))
    testSuite.addTest(test_nxprocess("test_align"))
    testSuite.addTest(test_nxprocess("test_expression"))
    #testSuite.addTest(test_nxprocess("test_resample"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

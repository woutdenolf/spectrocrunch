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
import itertools

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
        positioners.add_axis('x',range(shape[2])[::-1])
        yencres = 2
        xencres = 3
        ypossmax = 4
        xpossmax = 6

        posmap = np.arange(shape[1])[:,np.newaxis]+np.arange(shape[2])[np.newaxis,:]/(shape[2]-1.)*ypossmax
        yenc = np.stack([posmap*yencres]*shape[self.stackdim],axis=self.stackdim)
        posmap = np.arange(shape[2])[np.newaxis,::-1]+np.arange(shape[1])[:,np.newaxis]/(shape[1]-1.)*xpossmax
        xenc = np.stack([posmap*xencres]*shape[self.stackdim],axis=self.stackdim)

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
                           axis.factory(range(shape[2]+1)[::-1])
        elif method=='expression':
            info['expression'] = '{}/{arr_iodet}'
            info['skip'] = ['arr_iodet']
            info['select'] = process.results['counters']['arr_iodet']
        elif method=='resample':
            info['encoders'] = {'y':{'counter':'arr_samy','resolution':yencres},
                                'x':{'counter':'arr_samx','resolution':xencres}}
            info['shift'] = {'y':ypossmax,'x':xpossmax}
        
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
                    mi = np.min(data)*1.1
                    if mi==0:
                        mi = -1
                    data -= mi
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
    
    def _run_task(self,parameters,proc1):
        task = nxtask.newtask(previous=proc1,**parameters)
        task.run()
        proc2 = task.output
        self.assertTrue(task.done)
        task.run()
        proc3 = task.output
        self.assertEqual(proc2,proc3)
        task = nxtask.newtask(previous=proc1,**parameters)
        self.assertTrue(task.done)
        task.run()
        proc3 = task.output
        self.assertEqual(proc2,proc3)
        return proc2
        
    def test_crop(self):
        with self._nxprocess(method='crop') as proc1:
            proc1,info = proc1
            parameters = {'method':'crop','default':'Si-K','sliced':False,
                          'reference':'Al-K','nanval':np.nan}
            proc2 = self._run_task(parameters,proc1)
            
            parameters['sliced'] = True
            parameters['name'] = 'crop2'
            proc3 = self._run_task(parameters,proc1)
            self.assertNotEqual(proc2,proc3)
            
            grid1 = regulargrid.NXRegularGrid(proc1)
            grid2 = regulargrid.NXRegularGrid(proc2)
            grid3 = regulargrid.NXRegularGrid(proc3)
            self.assertEqual(set([sig.name for sig in grid1.signals]),
                             set([sig.name for sig in grid2.signals]))
            self.assertFalse(np.isnan(grid2.values).any())
            np.testing.assert_array_equal(grid2.values,grid3.values)
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
            proc2 = self._run_task(parameters,proc1)
            
            parameters['sliced'] = True
            parameters['name'] = 'replace2'
            proc3 = self._run_task(parameters,proc1)
            self.assertNotEqual(proc2,proc3)
            
            grid1 = regulargrid.NXRegularGrid(proc1)
            grid2 = regulargrid.NXRegularGrid(proc2)
            grid3 = regulargrid.NXRegularGrid(proc3)
            self.assertEqual(set([sig.name for sig in grid1.signals]),
                             set([sig.name for sig in grid2.signals]))
            np.testing.assert_array_equal(grid2.values,grid3.values)
            np.testing.assert_array_equal(grid1.values[info['index']],parameters['old'])
            np.testing.assert_array_equal(grid2.values[info['index']],parameters['new'])

            self.assertEqual(proc1.default.signal.name,parameters['default'])
    
    def test_minlog(self):
        with self._nxprocess(method='minlog') as proc1:
            proc1,info = proc1
            parameters = {'method':'minlog','sliced':False}
            proc2 = self._run_task(parameters,proc1)
            
            parameters['sliced'] = True
            parameters['name'] = 'minlog2'
            proc3 = self._run_task(parameters,proc1)
            self.assertNotEqual(proc2,proc3)
            
            grid1 = regulargrid.NXRegularGrid(proc1)
            grid2 = regulargrid.NXRegularGrid(proc2)
            grid3 = regulargrid.NXRegularGrid(proc3)
            
            self.assertEqual(set([sig.name for sig in grid1.signals]),
                             set([sig.name for sig in grid2.signals]))
            np.testing.assert_array_equal(grid2.values,grid3.values)
            np.testing.assert_array_equal(-np.log(grid1),grid3.values)
            
    def test_align(self):
        with self._nxprocess(method='align') as proc1:
            proc1,info = proc1
            parameters = {'method':'align','alignmethod':'max','reference':'Fe-K',
                          'refimageindex':0,'default':'Fe-K'}
            proc2 = self._run_task(parameters,proc1)
            
            grid2 = regulargrid.NXRegularGrid(proc2)
            axes = grid2.axes
            axes.pop(grid2.stackdim)
            for ax1,ax2 in zip(info['axes'],axes):
                self.assertEqual(ax1,ax2)
    
    def test_expression(self):
        with self._nxprocess(method='expression') as proc1:
            proc1,info = proc1
            skip = [{'method':'regex','pattern':name} for name in info['skip']]
            parameters = {'method':'expression','expression':info['expression'],'skip':skip,'sliced':False}
            proc2 = self._run_task(parameters,proc1)
            
            parameters['sliced'] = True
            parameters['name'] = 'expression2'
            proc3 = self._run_task(parameters,proc1)
            self.assertNotEqual(proc2,proc3)

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
            
            params = ((['y'],['x','y']),(True,False))
            
            for i,p in enumerate(itertools.product(*params),1):
                axes,crop = p
                encoders = {k:v for k,v in info['encoders'].items() if k in axes}
                parameters = {'name':'crop{}'.format(i),'method':'resample','encoders':encoders,'crop':crop}
                proc2 = self._run_task(parameters,proc1)

                grid1 = regulargrid.NXRegularGrid(proc1)
                grid2 = regulargrid.NXRegularGrid(proc2)
                
                # Check new axes position
                encoder_signals = {}
                offsets = proc2.results['encoder_offset'].read().tolist()
                offsets.insert(grid2.stackdim,0)
                for ax1,ax2,offset in zip(grid1.axes,grid2.axes,offsets):
                    name = ax1.name
                    if name in axes:
                        n = int(np.ceil(info['shift'][name]/2.))
                        if crop:
                            ax1 = axis.factory(ax1[n:-n])
                        else:
                            add = np.arange(1,n+1)*ax1.stepsize
                            addstart = ax1.start-add[::-1]
                            addend = ax1.end+add
                            x = addstart.magnitude.tolist()+\
                                ax1.magnitude.tolist()+\
                                addend.magnitude.tolist()
                            ax1 = axis.factory(units.Quantity(x,units=ax1.units))
                        resolution = units.Quantity(parameters['encoders'][name]['resolution'],units=1/ax2.units)
                        encoder_signals[name] = (ax2.values*resolution+offset).to('dimensionless').magnitude
                    self._check_axis(ax1,ax2)

                # Check encoder signals
                signals = grid2.signal_names
                for axname,encinfo in encoders.items():
                    i = signals.index(encinfo['counter'])
                    enc = grid2[i,...]
                    encnan = np.isnan(enc)
                    self.assertTrue(crop^encnan.any())
                    
                    # Expected encoder values
                    encvalues = encoder_signals[axname]
                    index = [np.newaxis]*enc.ndim
                    if axname=='y':
                        index[1] = slice(None)
                    else:
                        index[2] = slice(None)
                    encvalues = encvalues[index]
                    
                    # Handle nan's
                    if not crop:
                        m = np.ones(enc.shape)
                        m[encnan] = np.nan
                        encvalues = encvalues*m
                        encvalues[encnan] = 999
                        enc[encnan] = 999
                        
                    self.assertTrue(np.isclose(enc, encvalues).all())
                        
    def _check_axes(self,grid1,grid2):
        for ax1,ax2 in zip(grid1.axes,grid2.axes):
            self._check_axis(ax1,ax2)
    
    def _check_axis(self,ax1,ax2):
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
    testSuite.addTest(test_nxprocess("test_resample"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

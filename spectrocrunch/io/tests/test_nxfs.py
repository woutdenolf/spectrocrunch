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

import os
import unittest
from testfixtures import TempDirectory

from .. import fs
from .. import localfs
from .. import h5fs
from .. import nxfs
from ..utils import TemporaryFilename

class test_nxfs(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()
    
    def tearDown(self):
        self.dir.cleanup()
        
    def test_nxclasses(self):
        h5filename = os.path.join(self.dir.path,'test.h5')
        root = nxfs.Path('/',h5file=h5filename)
        
        nxroot = root.nxroot()
        
        self.assertRaises(ValueError,nxroot.nxentry,'entry0001',mode='r')
        entry1 = nxroot.nxentry('entry0001')
        self.assertEqual(entry1,root['entry0001'])
        self.assertEqual(entry1,nxroot['entry0001'])
        
        self.assertRaises(nxfs.NexusException,nxroot.nxsubentry,'subentrya')
        subentrya = entry1.nxsubentry('subentrya')
        self.assertEqual(entry1,subentrya.nxentry())
        
        self.assertRaises(nxfs.NexusException,nxroot.nxdata,'data1')
        data1 = subentrya.nxdata('data1')
        entry2 = data1.new_nxentry()
        self.assertEqual(entry2,root['entry0002'])

        self.assertRaises(ValueError,data1.nxinstrument,mode='r')
        instrument = data1.nxinstrument()
    
        self._check_nxdata(data1)
        self._check_process(entry2)

        #root.ls(recursive=True,stats=False,depth=3)
    
    def _check_process(self,entry):
        cfg1 = {'p1':10,'p2':[10,20]}
        process1,new = entry.nxprocess('fit',parameters=cfg1)
        self.assertTrue(new)
        _,new = entry.nxprocess('fit',parameters=cfg1)
        self.assertFalse(new)
        
        shape = (2,3)
        dtype = float
        signals = ['Fe-K','Si-K','Al-K','S-K','Ce-L']
        for detector in range(2):
            detector = process1.results['detector{:02d}'.format(detector)].mkdir()
            for name in signals:
                detector[name].mkfile(shape=shape,dtype=dtype)
        
        positioners = entry.positioners()
        positioners.add_axis('y',range(2),units='um',title='vertical')
        positioners.add_axis('x',range(3))
        
        for name in signals:
            process1.plotselect.add_signal(path=detector[name])
        process1.plotselect.set_axes('y','x')
        process1.plotselect.mark_default()
        
        self.assertEqual(process1.config.read(),cfg1)
        self.assertFalse(process1.previous_processes)
        
        process2,new = entry.nxprocess('align',previous=[process1])
        self.assertFalse(new)
        self.assertEqual(process2.config.read(),None)
        self.assertEqual(process2.previous_processes[0].linkdest(),process1)

        self.assertRaises(ValueError,entry.nxprocess,'align',parameters={'wrong':1},previous=[process1])
        
    def _check_nxdata(self,data1):
        y = 'y',range(2),{'units':'um','title':'vertical'}
        ywrong = 'y',[1,2],{'units':'um'}
        x = 'x',range(3),{}

        self.assertRaises(nxfs.NexusException, data1.set_axes, y,x)
        
        signals = ['Fe-K','Si-K','Al-K','S-K','Ce-L']
        for signal in signals:
            data1.add_signal(name=signal,dtype=int,shape=(len(y[1]),len(x[1])))
        self.assertRaises(fs.AlreadyExists,data1.add_signal,name=signals[-1],dtype=int,shape=(len(y[1])+1,len(x[1])))
        self._check_nxdata_signals(data1,signals[-1],signals[:-1])

        data1.remove_signal(signals[-1])
        self._check_nxdata_signals(data1,signals[-2],signals[:-2])
        
        data1.remove_signal(signals[0])
        self.assertEqual(data1.signal.name,signals[-2])
        self._check_nxdata_signals(data1,signals[-2],signals[1:-2])
        
        signal = signals[-2]
        signals = signals[1:-2]
        data1.default_signal(signals[0])
        self._check_nxdata_signals(data1,signals[0],signals[1:]+[signal])
        
        signal,signals = signals[0],signals[1:]+[signal]
        data1[signals[0]].mark_default()
        self._check_nxdata_signals(data1,signals[0],signals[1:]+[signal])

        data2 = data1.parent.nxdata('data2')
        for signal in data1.signals:
            data2.add_signal(path=signal)

        data1.set_axes(y,x)
        data1.set_axes(y,x)
        self.assertRaises(ValueError,data1.set_axes,ywrong,x)
        data2.set_axes('y','x')
        data3 = data1.parent['data3'].link('data2')
        self.assertEqual(data3.axes[0][1].units,'micrometer')
        
    def _check_nxdata_signals(self,data,signal,signals):
        self.assertEqual(data.signal.name,signal)
        for signal,name in zip(data.auxiliary_signals,signals):
            self.assertEqual(signal.name,name)
            
def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_nxfs("test_nxclasses"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
        

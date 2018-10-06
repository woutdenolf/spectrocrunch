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
from ..utils import TemporaryFilename

class test_fs(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()
    
    def tearDown(self):
        self.dir.cleanup()
        
    def test_local(self):
        path = os.path.join(self.dir.path,'local')
        root = localfs.Path(path)
        self.assertEqual(root,path)
        self.assertEqual(root['/a']['b']['c'].root,'/')
        self.assertEqual(root['a']['b']['c'].root,'/')
        self._check_path(root,root)
        print('')
        root.ls(recursive=True)
        
    def test_h5(self):
        h5filename = os.path.join(self.dir.path,'test.h5')
        root1 = h5fs.Path('/',h5file=h5filename)
        root2 = h5fs.Path(h5filename+':/')
        self.assertEqual(root1,h5filename+':/')
        self.assertEqual(root1,root2)
        h5filename = os.path.join(self.dir.path,'ext.h5')
        root2 = h5fs.Path(h5filename+':/')
        self._check_path(root1,root2,shape=(2,3),dtype=int)
        print('')
        root1.ls(recursive=True)
        
    def _check_path(self,root1,root2,**createparams):
        self.assertEqual(root1['a/b/c'].parent,root1['a/b'])
        self.assertEqual(root1['a/b/c/'].parent,root1['a/b'])
        self.assertEqual(root1['a/b/c.txt'].parent,root1['a/b'])
        
        # Create directory
        directory = root1['test_create']
        self.assertFalse(directory.exists)
        self.assertFalse(directory.isfile)
        self.assertFalse(directory.isdir)
        directory.mkdir()
        self.assertTrue(directory.exists)
        self.assertFalse(directory.isfile)
        self.assertTrue(directory.isdir)
        self.assertRaises(fs.AlreadyExists, directory.mkdir, force=False)
        directory.mkdir()
        self.assertTrue(directory.exists)
        self.assertFalse(directory.isfile)
        self.assertTrue(directory.isdir)

        # Create file
        file_atxt = directory['a.txt']
        self.assertFalse(file_atxt.exists)
        self.assertFalse(file_atxt.isfile)
        self.assertFalse(file_atxt.isdir)
        with self.assertRaises(fs.Missing):
            with file_atxt.open(mode='r+',**createparams):
                pass
        with file_atxt.open(mode='w',**createparams):
            pass
        
        self.assertTrue(file_atxt.exists)
        self.assertTrue(file_atxt.isfile)
        self.assertFalse(file_atxt.isdir)
        if createparams:
            with self.assertRaises(fs.AlreadyExists):
                with file_atxt.open(mode='x',**createparams):
                    pass
        with file_atxt.open(mode='r') as f:
            pass

        # Create Link
        dir_a = directory['a'].mkdir()
        lnk1 = directory['lnk1']
        lnk2 = directory['lnk2']
        lnk3 = directory['lnk3']
        lnk4 = directory['lnk4']

        lnk1.link(dir_a)
        lnk2.link(lnk1)
        lnk3.link(file_atxt)
        lnk4.link(lnk3)
        
        self.assertEqual(lnk1.linkdest(),dir_a)
        self.assertEqual(lnk2.linkdest(),lnk1)
        self.assertEqual(lnk3.linkdest(),file_atxt)
        self.assertEqual(lnk4.linkdest(),lnk3)
        
        self.assertEqual(lnk1.linkdest(follow=True),dir_a)
        self.assertEqual(lnk2.linkdest(follow=True),dir_a)
        self.assertEqual(lnk3.linkdest(follow=True),file_atxt)
        self.assertEqual(lnk4.linkdest(follow=True),file_atxt)

        lnk1["asub"].mkdir()
        with lnk1["asub.txt"].open(mode='w',**createparams):
            pass
        with directory['x.txt'].open(mode='w',**createparams) as f:
            pass
        lnk1["x.txt"].link(directory['x.txt'])
        
        # Create when wrong node is existing
        for path in [dir_a,lnk1,lnk2]:
            with self.assertRaises(fs.NotAFile):
                with path.open(mode='w',**createparams):
                    pass
                    
        for path in [file_atxt,lnk3,lnk4]:
            with self.assertRaises(fs.NotADirectory):
                path.mkdir()
        
        # Copy/Move
        lnk5 = lnk1.move(lnk1.parent['lnk5'])
        self.assertFalse(lnk1.exists)
        self.assertFalse(lnk2.linkdest().exists)
        
        file_btxt = file_atxt.copy(file_atxt.parent['b.txt'])
        file_ctxt = file_atxt.move(file_atxt.parent['c.txt'])
        self.assertFalse(lnk3.linkdest().exists)
        
        dir_b = dir_a.copy(dir_a.parent['b'])
        dir_c = dir_a.copy(dir_a.parent['c'],follow=True)
        dir_d = dir_c.move(dir_a.parent['d'])
        dir_e = dir_a.copy(dir_a.parent['e'],follow=True)
        self.assertFalse(dir_c.exists)
        self.assertTrue(dir_d.exists)
        self.assertTrue(dir_b.exists)
        self.assertTrue(dir_a.exists)
        self.assertTrue(dir_e.exists)
        self.assertNotEqual(dir_e['x.txt'],directory['x.txt'])
        self.assertEqual(dir_b['x.txt'].linkdest(follow=True),directory['x.txt'])
        
        lnkb = directory['lnkb']
        lnkd = directory['lnkd']
        lnkb.link(dir_b)
        lnkd.link(dir_d)
        
        lnkb.remove(recursive=False)
        lnkd.remove(recursive=True)
        self.assertFalse(lnkb.exists)
        self.assertFalse(lnkd.exists)
        self.assertTrue(dir_b.exists)
        self.assertFalse(dir_d.exists)

        self.assertFalse(lnk2.exists)
        self.assertFalse(lnk3.exists)
        self.assertFalse(lnk4.exists)
        self.assertTrue(lnk5.exists)
        
        with self.assertRaises(fs.DirectoryIsNotEmpty):
            dir_a.remove(recursive=False)
        
        # External
        data = root2['data3'].mkdir()
        with data['a.txt'].open(mode='w',**createparams):
            pass
        data['b.txt'].link(data['a.txt'],soft=True)
        data['c.txt'].link(data['a.txt'],soft=False)
        root2['data3'].copy(root2['data4'])
        
        root1['data1'].link(root2['data3'])
        root1['data2'].link(root2['data4'])
        
        root1['data2/a.txt'].copy(root1['data2a.txt'])
        root1['data2'].linkdest().copy(root1['data2b'])
        self.assertTrue(root1['data2b'].isdir)
        self.assertTrue(root1['data2a.txt'].isfile)
        
        root1['data2/a.txt'].move(root1['data2_a.txt'])
        root1['data1'].linkdest().move(root1['data1b'])
        self.assertFalse(root1['data1'].exists)
        self.assertFalse(root1['data2/a.txt'].exists)
        self.assertTrue(root1['data2b/a.txt'].exists)
        self.assertFalse(root2['data3'].exists)
        self.assertTrue(root1['data1b'].exists)
            
def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_fs("test_local"))
    testSuite.addTest(test_fs("test_h5"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
        

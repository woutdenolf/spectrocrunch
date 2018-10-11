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
from ..utils import TemporaryFilename

class test_localfs(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()
    
    def tearDown(self):
        self.dir.cleanup()
        
    def test_mode(self):
        # label   exist       not exist   read    write   position
        # r       -           error       -       error   start
        # r+      -           error       -       -       start
        # w       truncate    create      error   -       -
        # w+      truncate    create      -       -       -
        # x       error       create      error   -       -
        # x+      error       create      -       -       -
        # a       -           create      error   -       end
        # a+      -           create      -       -       end
        
        with TemporaryFilename(self.dir.path,suffix='.txt') as filename:
            # Doesn't exist
            self.assertRaises(fs.Missing, self._check_r, filename,'123')
            self._check_w(filename,'123')
            # Exist
            self._check_r(filename,'456')
            self._check_content(filename,'123')
            
        with TemporaryFilename(self.dir.path,suffix='.txt') as filename:
            # Doesn't exist
            self.assertRaises(fs.Missing, self._check_rp, filename,'123')
            self._check_w(filename,'123')
            # Exist
            self._check_rp(filename,'456')
            self._check_content(filename,'456')
            
        with TemporaryFilename(self.dir.path,suffix='.txt') as filename:
            # Doesn't exist
            self._check_w(filename,'123')
            self._check_content(filename,'123')
            # Exist
            self._check_w(filename,'456')
            self._check_content(filename,'456')
              
        with TemporaryFilename(self.dir.path,suffix='.txt') as filename:
            # Doesn't exist
            self._check_wp(filename,'123')
            self._check_content(filename,'123')
            # Exist
            self._check_wp(filename,'456')
            self._check_content(filename,'456')
        
        with TemporaryFilename(self.dir.path,suffix='.txt') as filename:
            # Doesn't exist
            self._check_a(filename,'123')
            self._check_content(filename,'123')
            # Exist
            self._check_a(filename,'456')
            self._check_content(filename,'123456')
        
        with TemporaryFilename(self.dir.path,suffix='.txt') as filename:
            # Doesn't exist
            self._check_ap(filename,'123')
            self._check_content(filename,'123')
            # Exist
            self._check_ap(filename,'456')
            self._check_content(filename,'123456')
            
        try:
            with TemporaryFilename(self.dir.path,suffix='.txt') as filename:
                # Doesn't exist
                self._check_x(filename,'123')
                self._check_content(filename,'123')
                # Exist
                self.assertRaises(IOError, self._check_x, filename, '456')
                
            with TemporaryFilename(self.dir.path,suffix='.txt') as filename:
                # Doesn't exist
                self._check_xp(filename,'123')
                self._check_content(filename,'123')
                # Exist
                self.assertRaises(IOError, self._check_xp, filename, '456')
        except ValueError:
            pass # Python 2 does not have x

    def _check_w(self,filename,word):
        with localfs.Path(filename,mode='w').open() as f:
            f.write(word)
            f.seek(0)
            self.assertRaises(IOError, f.read)

    def _check_wp(self,filename,word):
        with localfs.Path(filename,mode='w+').open() as f:
            f.write(word)
            f.seek(0)
            self.assertEqual(word,f.read())
    
    def _check_x(self,filename,word):
        with localfs.Path(filename,mode='x').open() as f:
            f.write(word)
            f.seek(0)
            self.assertRaises(IOError, f.read)
    
    def _check_xp(self,filename,word):
        with localfs.Path(filename,mode='x+').open() as f:
            f.write(word)
            f.seek(0)
            self.assertEqual(word,f.read())
    
    def _check_a(self,filename,word):
        with localfs.Path(filename,mode='a').open() as f:
            f.write(word)
            f.seek(0)
            self.assertRaises(IOError, f.read)

    def _check_ap(self,filename,word):
        with localfs.Path(filename,mode='a+').open() as f:
            f.write(word)
            n = len(word)
            f.seek(-n,1)
            self.assertEqual(word,f.read(n))
    
    def _check_r(self,filename,word):
        with localfs.Path(filename,mode='r').open() as f:
            self.assertRaises(IOError, f.write, word)
            f.read()
            
    def _check_rp(self,filename,word):
        with localfs.Path(filename,mode='r+').open() as f:
            f.write(word)
            f.seek(0)
            self.assertEqual(word,f.read(len(word)))
                  
    def _check_content(self,filename,word): 
        filename = str(filename)  
        if os.path.isfile(filename):
            with localfs.Path(filename,mode='r').open() as f:
                b = f.read()==word
        else:
            b = None==word
        self.assertTrue(b)

def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_localfs("test_mode"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
        

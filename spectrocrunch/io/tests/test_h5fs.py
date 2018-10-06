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

from .. import h5fs
from ..utils import TemporaryFilename

class test_h5fs(unittest.TestCase):

    def setUp(self):
        self.dir = TempDirectory()
    
    def tearDown(self):
        self.dir.cleanup()

    def test_mode(self):
        # label   exist       not exist   read    write
        # r       -           error       -       error
        # r+      -           error       -       -
        # w       truncate    create      -       -
        # x       error       create      -       -
        # a       -           create      -       -
        #
        # Comparision with OS: w, x and a are always with '+' (i.e. read as well as write)
        
        with TemporaryFilename(self.dir.path,suffix='.txt') as filename:
            # Doesn't exist
            self.assertRaises(IOError, self._check_r, filename,'123')
            self._check_w(filename,'123')
            # Exist
            self._check_r(filename,'456')
            self._check_content(filename,'123')
            
        with TemporaryFilename(self.dir.path,suffix='.txt') as filename:
            # Doesn't exist
            self.assertRaises(IOError, self._check_rp, filename,'123')
            self._check_w(filename,'123')
            # Exist
            self._check_rp(filename,'456')
            self._check_content(filename,'123456')

        with TemporaryFilename(self.dir.path,suffix='.txt') as filename:
            # Doesn't exist
            self._check_w(filename,'123')
            self._check_content(filename,'123')
            # Exist
            self._check_w(filename,'456')
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
            self._check_x(filename,'123')
            self._check_content(filename,'123')
            # Exist
            self.assertRaises(IOError, self._check_x, filename, '456')

    def _write(self,f,word):
        f.create_group(word)

    def _read(self,f,i=None):
        words = [k.split('/')[-1] for k in f]
        if i is not None:
            words = [words[-1]]
        return str(''.join(words))

    def _check_w(self,filename,word):
        with h5fs.h5File(filename,mode='w').open() as f:
            self._write(f,word)
            self.assertEqual(word,self._read(f))
    
    def _check_x(self,filename,word):
        with h5fs.h5File(filename,mode='x').open() as f:
            self._write(f,word)
            self.assertEqual(word,self._read(f))
    
    def _check_a(self,filename,word):
        with h5fs.h5File(filename,mode='a').open() as f:
            self._write(f,word)
            self.assertEqual(word,self._read(f,-1))
    
    def _check_r(self,filename,word):
        with h5fs.h5File(filename,mode='r').open() as f:
            self.assertRaises(ValueError, self._write, f, word)
            self._read(f)
            
    def _check_rp(self,filename,word):
        with h5fs.h5File(filename,mode='r+').open() as f:
            self._write(f,word)
                  
    def _check_content(self,filename,word): 
        filename = str(filename)  
        if os.path.isfile(filename):
            with h5fs.h5File(filename,mode='r').open() as f:
                b = self._read(f)==word
        else:
            b = None==word
        self.assertTrue(b)

def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_h5fs("test_mode"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
        

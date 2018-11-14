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
import numpy as np
from random import shuffle
import collections

from .. import hashing

class test_hashing(unittest.TestCase):

    def test_compare(self):
        eqfuncs = [hashing.hashequal,hashing.jhashequal,hashing.phashequal]
        a = {'a':1,'b':[1,{'c':2,'d':3}]}
        b = {'b':[1,{'c':2,'d':3}],'a':1}
        for func in eqfuncs:
            self.assertTrue(func(a,b))
    
        a = {'a':1,'b':[1,collections.OrderedDict([('c',2),('d',3)])]}
        b = {'b':[1,collections.OrderedDict([('c',2),('d',3)])],'a':1}
        for func in eqfuncs:
            self.assertTrue(func(a,b))
    
        b = {'b':[1,collections.OrderedDict([('d',3),('c',2)])],'a':1}
        for func in eqfuncs:
            self.assertFalse(func(a,b))
    
def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_hashing("test_compare"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

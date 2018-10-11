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
import os
import numpy as np

from ...utils import units
from .. import jsonpickle

def equal(a,b):
    if type(a)==np.ndarray:
        return np.array_equal(a,b)
    else:
        return a==b
        
class ExampleClass(object):

    def __init__(self,x,y,z=None):
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self,other):
        return equal(self.x,other.x) and equal(self.y,other.y) and equal(self.z,other.z)

    def __getstate__(self):
        return self.__dict__.copy()
    
    def __setstate__(self,state):
        self.__dict__.update(state)

class test_serialize(unittest.TestCase):

    def test_units(self):
        a = units.Quantity(5,"keV")
        b = 10.
        c = ExampleClass(10,[10,20],z=30)
        d = ExampleClass(np.array([1,2,3]),\
                ExampleClass(-10,-20,z=-30),\
                z=ExampleClass(-10,ExampleClass(-10,-20,z=-30),z=ExampleClass(-10,-20,z=-30)))
    
        data = {"a":a,"b":b,"c":c,"d":d}
        
        serialized = jsonpickle.encode(data)
        deserialized = jsonpickle.decode(serialized)
        
        self.assertEqual(data,deserialized)
        
def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_serialize("test_units"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

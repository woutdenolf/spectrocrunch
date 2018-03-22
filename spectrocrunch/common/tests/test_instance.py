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
import itertools

from .. import units
from .. import instance
from ...import ureg
from ...math import noisepropagation

import numpy as np

class test_instance(unittest.TestCase):

    def test_quantunits(self):
        for q,arr,expand in itertools.product([True,False],repeat=3):
            if arr:
                x = np.arange(1.,10)
            else:
                x = 10.

            if q:
                unit0 = ureg.millimeter
                if arr and expand:
                    y = np.vectorize(lambda a:units.Quantity(a,units=unit0),otypes=[object])(x)
                else:
                    y = units.Quantity(x,units=unit0)
                unit1 = ureg.meter
                munit = 1e-3
            else:
                unit0 = None
                unit1 = None
                munit = 1
                y = x

            # Test bubbling up and down
            if q:
                z = units.Quantity(x,units=unit0)
                np.testing.assert_array_equal(instance.asarray(y),instance.asarray(z))
                z = np.vectorize(lambda a:units.Quantity(a,units=unit0),otypes=[object])(x)
                np.testing.assert_array_equal(units.Quantity(y),units.Quantity(z))
                
            # Test magnitude
            a = x*munit
            b = units.magnitude(y,unit1)
            if arr:
                np.testing.assert_array_equal(a,b)
            else:
                self.assertEqual(a,b)

            # Test unit conversion
            if q:
                a = units.Quantity(x,units=unit0)
                b = units.quantity_like(y,a)
            if arr:
                np.testing.assert_array_equal(a,b)
            else:
                self.assertEqual(a,b)
    
    def _test_asarray(self,x,checktype=True):
        y,func = instance.asarrayf(x)
        self.assertTrue(isinstance(y,np.ndarray))
        y = func(y)
        if checktype:
            self.assertEqual(type(y),type(x))
             
    def test_asarray(self):
        a = np.int64(0)
        b = np.int64(1)
        nums = [a,np.array(a),[a],np.array([a]),[a,b],np.array([a,b])]
        checktypes = [True,True,False,True,False]
        for num,check in zip(nums,checktypes):
            x = num
            self._test_asarray(x,checktype=check)
            x = units.Quantity(num,units=ureg.millimeter)
            self._test_asarray(x)
            x = np.vectorize(lambda a:units.Quantity(a,units=ureg.millimeter),otypes=[object])(x)
            self._test_asarray(x,checktype=False)
            x = noisepropagation.randomvariable(num,num)
            self._test_asarray(x,checktype=check)
        
        objs = ["abc",np.array("abc"),np.array(["abc"]),["abc"],("abc",1,2)]
        checktypes = [True,True,True,False,False]
        for o,check in zip(nums,checktypes):
            self._test_asarray(x,checktype=check)

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_instance("test_asarray"))
    #testSuite.addTest(test_instance("test_quantunits"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
        

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

from .. import axis
from ...utils import units

class test_axis(unittest.TestCase):

    def test_quantitative(self):
        a = units.Quantity(10,units='um')
        i = units.Quantity(13,units='um')
        j = units.Quantity(17,units='um')
        b = units.Quantity(20,units='um')
        precision = units.Quantity(1,units='nm')
        n = 10
        
        ax1 = axis.AxisRegular(a,b,n,title='x',precision=precision)
        self.assertEqual(ax1,axis.Axis(ax1.values))
        
        ax2 = axis.AxisRegularInc(a.to('mm'),ax1.stepsize,n,title='y',precision=precision)
        self.assertEqual(ax2,axis.Axis(ax2.values))
        self.assertEqual(ax1,ax2)

        ax3 = axis.AxisNumber(a,title='z',precision=precision)
        self.assertEqual(ax3,axis.Axis(ax3.values))
        self.assertEqual(ax3,axis.AxisRegular(a,a,1))
        
        ax4 = axis.Axis(ax2.values,precision=ax2.precision).simplify()
        self.assertTrue(isinstance(ax4,axis.AxisRegular))
        
        o = units.Quantity(3,units='um')
        ax5 = axis.AxisSegments(a,n,i,n,j,n,b)
        ax6 = axis.AxisSegments(a-o,n,i,n,j,n,b)
        ax5.limits = a-o,i,j,b
        self.assertEqual(ax5,ax6)
        
        
def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_axis("test_quantitative"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

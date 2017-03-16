# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
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
from spectrocrunch.align.tests import test_all as test_align
from spectrocrunch.common.tests import test_all as test_common
from spectrocrunch.io.tests import test_all as test_io
from spectrocrunch.materials.tests import test_all as test_materials
from spectrocrunch.math.tests import test_all as test_math
from spectrocrunch.process.tests import test_all as test_process
from spectrocrunch.visualization.tests import test_all as test_visualization
from spectrocrunch.xrf.tests import test_all as test_xrf
from spectrocrunch.h5stacks.tests import test_all as test_h5stacks

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_align.test_suite_all())
    testSuite.addTest(test_common.test_suite_all())
    testSuite.addTest(test_io.test_suite_all())
    testSuite.addTest(test_materials.test_suite_all())
    testSuite.addTest(test_math.test_suite_all())
    testSuite.addTest(test_process.test_suite_all())
    testSuite.addTest(test_visualization.test_suite_all())
    testSuite.addTest(test_xrf.test_suite_all())
    testSuite.addTest(test_h5stacks.test_suite_all())
    return testSuite

if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

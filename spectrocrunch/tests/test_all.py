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

from ..patch.tests import test_all as test_patch
from ..align.tests import test_all as test_align
from ..utils.tests import test_all as test_utils
from ..detectors.tests import test_all as test_detectors
from ..fullfield.tests import test_all as test_fullfield
from ..geometries.tests import test_all as test_geometries
from ..h5stacks.tests import test_all as test_h5stacks
from ..instruments.tests import test_all as test_instruments
from ..io.tests import test_all as test_io
from ..materials.tests import test_all as test_materials
from ..math.tests import test_all as test_math
from ..process.tests import test_all as test_process
from ..optics.tests import test_all as test_optics
from ..pipelines.tests import test_all as test_pipeline
from ..simulation.tests import test_all as test_simulation
from ..sources.tests import test_all as test_sources
from ..visualization.tests import test_all as test_visualization
from ..xrf.tests import test_all as test_xrf


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_patch.test_suite())
    testSuite.addTest(test_utils.test_suite())
    testSuite.addTest(test_math.test_suite())
    testSuite.addTest(test_io.test_suite())
    testSuite.addTest(test_instruments.test_suite())
    testSuite.addTest(test_sources.test_suite())
    testSuite.addTest(test_geometries.test_suite())
    testSuite.addTest(test_detectors.test_suite())
    testSuite.addTest(test_process.test_suite())
    testSuite.addTest(test_optics.test_suite())
    testSuite.addTest(test_materials.test_suite())
    testSuite.addTest(test_simulation.test_suite())
    testSuite.addTest(test_align.test_suite())
    testSuite.addTest(test_fullfield.test_suite())
    testSuite.addTest(test_xrf.test_suite())
    testSuite.addTest(test_h5stacks.test_suite())
    testSuite.addTest(test_pipeline.test_suite())
    testSuite.addTest(test_visualization.test_suite())
    return testSuite

if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

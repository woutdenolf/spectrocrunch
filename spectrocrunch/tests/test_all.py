# -*- coding: utf-8 -*-

import unittest

from ..testutils.tests import test_all as test_testutils
from ..patch.tests import test_all as test_patch
from ..align.tests import test_all as test_align
from ..utils.tests import test_all as test_utils
from ..detectors.tests import test_all as test_detectors
from ..fullfield.tests import test_all as test_fullfield
from ..geometries.tests import test_all as test_geometries
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
    testSuite.addTest(test_testutils.test_suite())
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
    testSuite.addTest(test_pipeline.test_suite())
    testSuite.addTest(test_visualization.test_suite())
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

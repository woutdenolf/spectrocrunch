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

from ..import base
from ...patch import jsonpickle


class test_base(unittest.TestCase):

    def test_serialize(self):
        g1 = base.Base()
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)

        g1 = base.SolidAngle()
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)
        g1 = base.SolidAngle(solidangle=0.8)
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)

        g1 = base.FlatSample()
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)
        g1 = base.FlatSample(anglein=40.1, angleout=-50., azimuth=3.)
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)

        g1 = base.Centric()
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)
        g1 = base.Centric(distance=0.1, anglein=40.1,
                          angleout=-50., azimuth=3.)
        g2 = jsonpickle.loads(jsonpickle.dumps(g1))
        self.assertEqual(g1, g2)


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_base("test_serialize"))
    return testSuite


if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

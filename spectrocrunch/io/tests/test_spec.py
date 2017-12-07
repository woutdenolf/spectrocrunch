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

from .. import spec

class test_spec(unittest.TestCase):

    def test_cmd_parser(self):
        cmds = [("zapimage sampy 0 1 10 sampz 2 3 11 100",\
                {"time":lambda x: x.to("ms").magnitude==100}),\
                ("zapimage sampy 0 1 10 100 sampz 2 3 11",\
                {"time":lambda x: x.to("ms").magnitude==100}),\
                ("puzzle sampy 0 1 10 sampz 2 3 11 100",\
                {"time":lambda x: x.to("ms").magnitude==100}),\
                ("mesh sampy 0 1 10 sampz 2 3 11 0.1",\
                {"time":lambda x: x.to("ms").magnitude==100}),\
                ("zapline sampy 0 1 10 100",\
                {"time":lambda x: x.to("ms").magnitude==100}),\
                ("ascan sampy 0 1 10 0.1",\
                {"time":lambda x: x.to("ms").magnitude==100}),\
                ("zapenergy SUM 10 100",\
                {"time":lambda x: x.to("ms").magnitude==100}),\
                ("zapenergy SUM2 10 100",\
                {"time":lambda x: x.to("ms").magnitude==100}),\
                ("invalid",{"name":lambda x: x=="unknown"})]

        p = spec.cmd_parser()
        for cmd,checks in cmds:
            r = p.parse(cmd)
            for k,func in checks.items():
                self.assertTrue(func(r[k]))

def test_suite_all():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_spec("test_cmd_parser"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite_all()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)
        

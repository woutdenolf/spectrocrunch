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
import os
import sys

from .. import signalhandling

class test_signalhandling(unittest.TestCase):

    def mysetup(self):
        self.state = 1

    def myteardown(self, exc_type, exc_value, exc_traceback):
        if exc_type:
            self.state = 0
        return 1 # signal is not propagated

    def myteardown_propagate(self, exc_type, exc_value, exc_traceback):
        pass

    def _check_signal(self,sendsignal):
        with signalhandling.DelaySignalsContext(setup=self.mysetup,teardown=self.myteardown):
            self.state += 1
            sendsignal()
            self.state += 1
        self.assertTrue(self.state==0)

    def test_noerror(self):
        with signalhandling.DelaySignalsContext(setup=self.mysetup,teardown=self.myteardown):
            self.state += 1
        self.assertTrue(self.state==2)

    def test_error(self):
        with signalhandling.DelaySignalsContext(setup=self.mysetup,teardown=self.myteardown):
            self.state += 1
            raise RuntimeError()
            self.state += 1
        self.assertTrue(self.state==0)

        with self.assertRaises(RuntimeError):
            with signalhandling.DelaySignalsContext(setup=self.mysetup,teardown=self.myteardown_propagate):
                self.state += 1
                raise RuntimeError()
                self.state += 1
        self.assertTrue(self.state==2)

    def test_sigterm(self):
        self._check_signal(lambda : os.kill(os.getpid(), signal.SIGTERM))
    
    def test_sigint(self):
        self._check_signal(lambda : os.kill(os.getpid(), signal.SIGINT))

    def test_sigexit(self):
        self._check_signal(lambda : sys.exit(0))
        self._check_signal(lambda : exit(0))

    def test_sigexit(self):
        self._check_signal(lambda : sys.exit(0))
        self._check_signal(lambda : exit(0))


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_signalhandling("test_noerror"))
    testSuite.addTest(test_signalhandling("test_error"))
    testSuite.addTest(test_signalhandling("test_sigterm"))
    testSuite.addTest(test_signalhandling("test_sigint"))
    testSuite.addTest(test_signalhandling("test_sigexit"))
    return testSuite
    
if __name__ == '__main__':
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)


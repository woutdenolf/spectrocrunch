# -*- coding: utf-8 -*-

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
        return 1  # signal is not propagated

    def myteardown_propagate(self, exc_type, exc_value, exc_traceback):
        pass

    def _check_signal(self, sendsignal):
        with signalhandling.HandleTermination(setup=self.mysetup, teardown=self.myteardown):
            self.state += 1
            sendsignal()
            self.state += 1
        self.assertTrue(self.state == 0)

    def test_noerror(self):
        with signalhandling.HandleTermination(setup=self.mysetup, teardown=self.myteardown):
            self.state += 1
        self.assertTrue(self.state == 2)

    def test_error(self):
        with signalhandling.HandleTermination(setup=self.mysetup, teardown=self.myteardown):
            self.state += 1
            raise RuntimeError()
            self.state += 1
        self.assertTrue(self.state == 0)

        with self.assertRaises(RuntimeError):
            with signalhandling.HandleTermination(setup=self.mysetup, teardown=self.myteardown_propagate):
                self.state += 1
                raise RuntimeError()
                self.state += 1
        self.assertTrue(self.state == 2)

    def test_sigterm(self):
        self._check_signal(lambda: os.kill(os.getpid(), signal.SIGTERM))

    def test_sigint(self):
        self._check_signal(lambda: os.kill(os.getpid(), signal.SIGINT))

    def test_sigexit(self):
        self._check_signal(lambda: sys.exit(0))
        self._check_signal(lambda: exit(0))

    def test_sigexit(self):
        self._check_signal(lambda: sys.exit(0))
        self._check_signal(lambda: exit(0))


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

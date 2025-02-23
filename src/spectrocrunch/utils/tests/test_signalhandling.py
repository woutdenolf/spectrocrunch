import unittest
import os
import sys
import signal

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
        with signalhandling.HandleTermination(
            setup=self.mysetup, teardown=self.myteardown
        ):
            self.state += 1
            sendsignal()
            # self.state += 1  # this makes test_sigterm fail
        self.assertTrue(self.state == 0)

    def test_noerror(self):
        with signalhandling.HandleTermination(
            setup=self.mysetup, teardown=self.myteardown
        ):
            self.state += 1
        self.assertTrue(self.state == 2)

    def test_error(self):
        with signalhandling.HandleTermination(
            setup=self.mysetup, teardown=self.myteardown
        ):
            self.state += 1
            raise RuntimeError()
            self.state += 1
        self.assertTrue(self.state == 0)

        with self.assertRaises(RuntimeError):
            with signalhandling.HandleTermination(
                setup=self.mysetup, teardown=self.myteardown_propagate
            ):
                self.state += 1
                raise RuntimeError()
                self.state += 1
        self.assertTrue(self.state == 2)

    @unittest.skipIf(sys.platform == "win32", "Skipping on Windows")
    def test_sigterm(self):
        self._check_signal(lambda: os.kill(os.getpid(), signal.SIGTERM))

    @unittest.skipIf(sys.platform == "win32", "Skipping on Windows")
    def test_sigint(self):
        with self.assertRaises(KeyboardInterrupt):
            self._check_signal(lambda: os.kill(os.getpid(), signal.SIGINT))

    @unittest.skipIf(sys.platform == "win32", "Skipping on Windows")
    def test_sigexit(self):
        self._check_signal(lambda: sys.exit(0))
        self._check_signal(lambda: exit(0))

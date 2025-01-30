# -*- coding: utf-8 -*-

import os
import unittest
from testfixtures import TempDirectory

from .. import fs
from .. import localfs
from .. import h5fs
from .. import nxfs
from ..utils import TemporaryFilename


class test_nxfs(unittest.TestCase):
    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()

    def test_nxclasses(self):
        h5filename = os.path.join(self.dir.path, "test.h5")
        root = nxfs.Path("/", h5file=h5filename)

        nxroot = root.nxroot()

        self.assertRaises(ValueError, nxroot.nxentry, "entry0001", mode="r")
        entry1 = nxroot.nxentry("entry0001")
        self.assertEqual(entry1, root["entry0001"])
        self.assertEqual(entry1, nxroot["entry0001"])

        self.assertRaises(nxfs.NexusException, nxroot.nxsubentry, "subentrya")
        subentrya = entry1.nxsubentry("subentrya")
        self.assertEqual(entry1, subentrya.nxentry())

        self.assertRaises(nxfs.NexusException, nxroot.nxdata, "data1")
        data1 = subentrya.nxdata("data1")
        entry2 = data1.new_nxentry()
        self.assertEqual(entry2, root["entry0002"])

        self.assertRaises(ValueError, data1.nxinstrument, mode="r")
        instrument = data1.nxinstrument()

        self._check_nxdata(data1)
        self._check_process(entry2)

        # root.ls(recursive=True,stats=False,depth=3)

    def _check_process(self, entry):
        cfg1 = {"p1": 10, "p2": [10.0, 20.0], "p3": {"a": 1}, "p4": "test"}

        process1 = entry.find_nxprocess(parameters=cfg1)
        self.assertIs(process1, None)
        process1 = entry.nxprocess("fit", parameters=cfg1)
        self.assertTrue(process1.exists)
        process1 = entry.find_nxprocess(parameters=cfg1)
        self.assertIsNot(process1, None)
        process1 = entry.nxprocess("fit", parameters=cfg1)
        self.assertTrue(process1.exists)

        shape = (2, 3)
        dtype = float
        signals = ["Fe-K", "Si-K", "Al-K", "S-K", "Ce-L"]
        for detector in range(2):
            detector = process1.results["detector{:02d}".format(detector)].mkdir()
            for name in signals:
                detector[name].mkfile(shape=shape, dtype=dtype)

        positioners = process1.positioners()
        positioners.add_axis("y", range(2), units="um", title="vertical")
        positioners.add_axis("x", range(3))

        for name in signals:
            process1.plotselect.add_signal(path=detector[name])
        process1.plotselect.set_axes("y", "x")
        process1.plotselect.mark_default()

        self.assertEqual(process1.config.read(), cfg1)
        self.assertFalse([dep for dep in process1.dependencies])

        process2 = entry.find_nxprocess(dependencies=[process1])
        self.assertFalse(process2)
        process2 = entry.nxprocess("process", dependencies=[process1])
        self.assertTrue(process2.exists)
        self.assertEqual(process2.config.read(), None)
        self.assertEqual(next(iter(process2.dependencies)).linkdest(), process1)

        process2b = entry.find_nxprocess(dependencies=[process1])
        self.assertEqual(process2, process2b)

        process2b = entry.find_nxprocess(
            dependencies=[process1], parameters={"wrong": 1}
        )
        self.assertIs(process2b, None)
        process2b = entry.nxprocess(
            "process", dependencies=[process1], parameters={"wrong": 1}
        )
        self.assertTrue(process2b.exists)
        self.assertNotEqual(process2, process2b)

    def _check_nxdata(self, data1):
        y = "y", range(2), {"units": "um", "title": "vertical"}
        ywrong = "y", [1, 2], {"units": "um"}
        x = "x", range(3), {}

        signals = ["Fe-K", "Si-K", "Al-K", "S-K", "Ce-L"]
        for signal in signals:
            data1.add_signal(name=signal, dtype=int, shape=(len(y[1]), len(x[1])))
        self.assertRaises(
            nxfs.NexusFormatException,
            data1.add_signal,
            name=signals[-1],
            dtype=int,
            shape=(len(y[1]) + 1, len(x[1])),
        )
        self.assertRaises(
            fs.AlreadyExists,
            data1.add_signal,
            name=signals[-1],
            dtype=int,
            shape=(len(y[1]), len(x[1])),
        )
        self._check_nxdata_signals(data1, signals[-1], signals[:-1])

        data1.remove_signal(signals[-1])
        self._check_nxdata_signals(data1, signals[-2], signals[:-2])

        data1.remove_signal(signals[0])
        self.assertEqual(data1.signal.name, signals[-2])
        self._check_nxdata_signals(data1, signals[-2], signals[1:-2])

        signal = signals[-2]
        signals = signals[1:-2]
        data1.default_signal(signals[0])
        self._check_nxdata_signals(data1, signals[0], [signal] + signals[1:])

        signal, signals = signals[0], signals[1:] + [signal]
        data1[signals[0]].mark_default()
        self._check_nxdata_signals(data1, signals[0], [signal] + signals[1:])

        data2 = data1.parent.nxdata("data2")
        for signal in data1.signals:
            data2.add_signal(path=signal)

        data1.set_axes(y, x)
        data1.set_axes(y, x)
        self.assertRaises(ValueError, data1.set_axes, ywrong, x)
        data2.set_axes("y", "x")
        data3 = data1.parent["data3"].link("data2")
        self.assertEqual(data3.axes[0][1].units, "micrometer")

    def _check_nxdata_signals(self, data, signal, signals):
        self.assertEqual(data.signal.name, signal)
        for signal, name in zip(data.auxiliary_signals, signals):
            self.assertEqual(signal.name, name)


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_nxfs("test_nxclasses"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

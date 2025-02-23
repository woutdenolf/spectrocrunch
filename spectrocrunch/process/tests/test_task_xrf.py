import unittest
import os

from ...io import nxfs
from ...patch import xraylib
from .xrfmap import XrfMapGenerator
from .test_task import test_task


class test_task_xrf(test_task):
    @unittest.skipIf(xraylib.XRayInit is None, "xraylib not installed")
    def test_nxqxrf(self):
        fixed = {
            "plot": False,
            "dark": False,
            "time": 1,
            "gaindiodeI0": 1e8,
            "gaindiodeIt": 1e7,
        }
        variable = [
            {"I0_counts": 300, "It_counts": 30, "dark": True},
            {"I0_counts": 400000, "It_counts": 100000, "energy": 7},
            {"I0_counts": 200000, "It_counts": 100000, "energy": 7.1},
        ]
        h5filename = os.path.join(self.dir.path, "test.h5")
        outputparent = nxfs.Path("/", h5file=h5filename)["entry"]
        parameters = {
            "method": "xrfgeometry",
            "outputparent": outputparent,
            "geometry": "sxm1",
            "init": {},
            "fixed": fixed,
            "variable": variable,
        }
        proc2 = self._run_task(parameters, None)
        parameters["geometry"] = "SXM1"
        proc3 = self._run_task(parameters, None)
        self._check_reproc(proc2, proc3)

    @unittest.skipIf(xraylib.XRayInit is None, "xraylib not installed")
    def test_nxxiaedf(self):
        xrfmap = XrfMapGenerator(nmaps=1)
        xrfmap.generate(self.dir.path, "test")

        # Needs a quantification run
        h5filename = os.path.join(self.dir.path, "geometries.h5")
        seldetectors = [(det,) for det in range(xrfmap.ndet)]
        parameters = xrfmap.taskparams_geometry(seldetectors)
        parameters["method"] = "xrfgeometry"
        parameters["outputparent"] = h5filename + "::/xrf"
        proc1 = self._run_task(parameters, None)

        # Output is an NXentry (not an NXprocess)
        h5filename = os.path.join(self.dir.path, "test.h5")
        parameters = {
            "method": "xiaedftonx",
            "name": "first",
            "outputparent": h5filename,
            "path": xrfmap.path,
            "radix": xrfmap.radix,
            "number": xrfmap.scannumbers[0],
            "instrument": xrfmap.instrument.lower(),
        }
        parameters.update(xrfmap.taskparams_pymca(seldetectors))
        entry2 = self._run_task(parameters, proc1, outputnxprocess=False)
        parameters["name"] = "second"
        entry3 = self._run_task(parameters, proc1, outputnxprocess=False)
        appli2 = entry2.root.find_application(entryname="first", definition="NXxrf")
        appli3 = entry3.root.find_application(entryname="second", definition="NXxrf")
        self.assertEqual(entry2, appli2.parent)
        self.assertEqual(entry3, appli3.parent)
        files2 = {path.name for path in appli2.listdir()}
        files3 = {path.name for path in appli3.listdir()}
        files = {
            "i0",
            "i0_to_flux_factor",
            "i0_to_flux_offset",
            "it",
            "it_to_flux_factor",
            "it_to_flux_offset",
            "start_time",
            "end_time",
            "definition",
        }
        files |= {"mca{:02d}".format(i) for i in range(xrfmap.ndet)}
        self.assertEqual(files, files2)
        self.assertEqual(files, files3)


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_task_xrf("test_nxqxrf"))
    testSuite.addTest(test_task_xrf("test_nxxiaedf"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

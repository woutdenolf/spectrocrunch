# -*- coding: utf-8 -*-

import unittest
from testfixtures import TempDirectory
import contextlib
import os
import numpy as np
import itertools
import logging
import sys

from .. import id21_ffxas
from ..run import run_sequential
from ...io.edf import saveedf
from ...process.h5regulargrid import NXSignalRegularGrid

logger = logging.getLogger(__name__)


class test_ffxas(unittest.TestCase):
    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()

    def _data_generate(self, path, radix):
        n1, n2 = 55, 65
        n = 3

        xv = range(int(n2 / 3), int(n2 / 3) + n)
        yv = range(int(n2 / 3), int(n2 / 3) + n)
        energy = np.linspace(7, 7.3, n)  # keV
        intensity = np.linspace(200, 100, n)  # DU/sec
        transmission = np.linspace(0.9, 0.3, n)
        nbdata = 8
        nbflat = 2
        tdata = 0.4  # sec
        tflat = 0.1  # sec

        darkoff = 50  # DU
        darkgain = 33  # DU/sec
        nbdark = 10

        darkdata = np.full((n1, n2), (darkoff + darkgain * tdata) * nbdark)
        hdark = {"energy": energy.max(), "exposure_time": tdata, "nb_frames": nbdark}
        filename = os.path.join(path, "{}_dark_{}_0000.edf".format(radix, tdata))
        saveedf(filename, darkdata, hdark)

        if tflat != tdata:
            darkflat = np.full((n1, n2), (darkoff + darkgain * tflat) * nbdark)
            hdark = {
                "energy": energy.max(),
                "exposure_time": tflat,
                "nb_frames": nbdark,
            }
            filename = os.path.join(path, "{}_dark_{}_0000.edf".format(radix, tflat))
            saveedf(filename, darkflat, hdark)

        for c, (x, y, e, i, t) in enumerate(
            zip(xv, yv, energy, intensity, transmission)
        ):
            data = np.full((n1, n2), (darkoff + (i + darkgain) * tdata) * nbdata)
            data[y, x] = (darkoff + (t * i + darkgain) * tdata) * nbdata

            hdata = {"energy": e, "exposure_time": tdata, "nb_frames": nbdata}
            filename = os.path.join(
                path, "{}_data_0000_{:04d}_0000.edf".format(radix, c)
            )
            saveedf(filename, data, hdata)

            flat1 = np.full(
                (n1, n2), (darkoff + (i * 1.05 + darkgain) * tflat) * nbflat
            )
            flat2 = np.full(
                (n1, n2), (darkoff + (i * 0.95 + darkgain) * tflat) * nbflat
            )

            hflat = {"energy": e, "exposure_time": tflat, "nb_frames": nbflat}
            filename = os.path.join(
                path, "{}_ref_0000_{:04d}_{:04d}.edf".format(radix, c, 0)
            )
            saveedf(filename, flat1, hflat)
            filename = os.path.join(
                path, "{}_ref_0000_{:04d}_{:04d}.edf".format(radix, c, 1)
            )
            saveedf(filename, flat2, hflat)

        self._datainfo = {
            "sourcepaths": path,
            "radix": radix,
            "energy": energy,
            "intensity": intensity,
            "transmission": transmission,
            "xv": xv,
            "yv": yv,
            "n1": n1,
            "n2": n2,
            "n": n,
        }

    @contextlib.contextmanager
    def _destpath_context(self):
        destpath = TempDirectory()
        yield destpath
        destpath.cleanup()

    def _process(self, crop, roialign, normalizeonload, stackdim):
        parameters = {"nxentry": "entry_0001"}

        # Raw data
        parameters["sourcepaths"] = self._datainfo["sourcepaths"]
        radix = self._datainfo["radix"]
        parameters["radix"] = radix
        parameters["rebin"] = (1, 1)
        parameters["stackdim"] = stackdim

        # Normalization
        parameters["normalize"] = True
        parameters["normalizeonload"] = normalizeonload

        # Alignment
        parameters["alignmethod"] = "max"
        parameters["refimageindex"] = -1
        parameters["plot"] = False
        parameters["crop"] = crop
        parameters["roiraw"] = None
        parameters["roialign"] = roialign
        parameters["roiresult"] = None

        with self._destpath_context() as destpath:
            parameters["nxentry"] = (
                os.path.join(destpath.path, radix + ".h5") + "::/" + radix
            )
            for repeat in range(2):
                logger.debug(
                    "Process parameters\n"
                    + "\n".join(["{}:{}".format(k, v) for k, v in parameters.items()])
                )

                tasks = id21_ffxas.tasks(**parameters)
                if repeat:
                    for task in tasks:
                        self.assertTrue(task.done)
                else:
                    for task in tasks:
                        self.assertFalse(task.done)
                    run_sequential(tasks)
                    for task in tasks:
                        self.assertTrue(task.done)
                    self._checkresult(parameters, tasks)

    def _checkresult(self, parameters, tasks):
        params = self._datainfo
        imgshape = (params["n1"], params["n2"])
        shape = [params["n1"], params["n2"]]
        shape.insert(parameters["stackdim"], params["n"])
        shape = tuple(shape)

        # Check nxprocess groups
        groups = ["process:fullfield.1"]
        if parameters["normalize"] and not parameters["normalizeonload"]:
            groups.append("process:normalize.1")
        if parameters["alignmethod"]:
            groups.append("process:align.1")
        if parameters["roiresult"]:
            groups.append("process:roi.1")
        entry = tasks[-1].output.nxentry()
        self.assertEqual(
            set(groups), set([g.name for g in entry.iter_is_nxclass("NXprocess")])
        )

        # Check axes
        nxprocess = entry["process:fullfield.1"]
        positioners = nxprocess.positioners()
        np.testing.assert_allclose(params["energy"], positioners.get_axis("energy"))
        np.testing.assert_array_equal(range(params["n1"]), positioners.get_axis("row"))
        np.testing.assert_array_equal(range(params["n2"]), positioners.get_axis("col"))

        # Check transmission
        if parameters["normalize"]:
            if parameters["normalizeonload"]:
                nxprocess = entry["process:fullfield.1"]
            else:
                nxprocess = entry["process:normalize.1"]
            fdata = NXSignalRegularGrid(nxprocess.results["detector0"].signal)
            self.assertEqual(shape, fdata.shape)
            index = [slice(None)] * fdata.ndim
            for i, (transmission, x, y) in enumerate(
                zip(params["transmission"], params["xv"], params["yv"])
            ):
                index = [y, x]
                index.insert(parameters["stackdim"], params["energy"][i])
                index = fdata.locate(*index)
                transmission_calc = np.exp(-fdata[index])
                np.testing.assert_allclose(transmission_calc, transmission, rtol=1e-5)

        # Check aligned results
        if parameters["alignmethod"]:
            nxprocess = entry["process:align.1"]
            fdata = NXSignalRegularGrid(nxprocess.results["detector0"].signal)

            i = parameters["refimageindex"]
            x = params["xv"][i]
            y = params["yv"][i]
            for i, transmission in enumerate(params["transmission"]):
                index = [y, x]
                index.insert(parameters["stackdim"], params["energy"][i])
                index = fdata.locate(*index)
                transmission_calc = np.exp(-fdata[index])
                np.testing.assert_allclose(transmission_calc, transmission, rtol=1e-5)

    def test_process(self):
        sourcepaths = self.dir.path
        radix = "ff"
        self._data_generate(sourcepaths, radix)
        parameters = [(True, False), (None, ((3, -3), (4, -4))), (True, False), (0,)]
        if hasattr(self, "subTest"):
            for i, combination in enumerate(itertools.product(*parameters)):
                with self.subTest(i=i):
                    self._process(*combination)
                    sys.stdout.write(".")
                    sys.stdout.flush()
        else:
            for i, combination in enumerate(itertools.product(*parameters)):
                self._process(*combination)
                sys.stdout.write(".")
                sys.stdout.flush()


def test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_ffxas("test_process"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

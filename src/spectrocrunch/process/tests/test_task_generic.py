import unittest
import os
from contextlib import contextmanager
import itertools
import threading
import numpy as np

from .test_task import test_task
from .. import h5regulargrid
from ...io import nxfs
from ...utils.tests import genindexing
from .. import utils
from .. import axis
from ...utils import units
from ...io import target


class test_task_generic(test_task):
    @contextmanager
    def _nxprocess(self, method=None):
        h5filename = os.path.join(self.dir.path, "test.h5")
        root = nxfs.Path("/", h5file=h5filename).nxroot()
        entry = root.new_nxentry()
        parameters = {"name": "fromraw", "a": 1, "b": 2}
        name = target.prepare(parameters["name"])
        nxprocess = entry.nxprocess(name, parameters=parameters, dependencies=None)
        info = {}

        shape = (2, 10, 13)
        self.stackdim = 0
        nstack, nhor, nvert = shape
        nstack = nstack
        z = "z", range(nvert), {"units": "um", "title": "vertical"}
        y = "y", range(nhor), {"units": "um", "title": "horizontal"}
        x = "x", range(nstack), None
        yencres = 2
        zencres = 3
        ypossmax = 4
        zpossmax = 6

        posmap = (
            np.arange(nhor)[np.newaxis, :]
            + np.arange(nvert)[:, np.newaxis] / (nvert - 1.0) * ypossmax
        )
        yenc = np.stack([posmap.T * yencres] * nstack, axis=self.stackdim)
        posmap = (
            np.arange(nvert)[:, np.newaxis]
            + np.arange(nhor)[np.newaxis, :] / (nhor - 1.0) * zpossmax
        )
        zenc = np.stack([posmap.T * zencres] * nstack, axis=self.stackdim)

        dtype = np.float32
        signals = ["Fe-K", "Si-K", "Al-K", "S-K", "Ce-L"]
        counters = ["arr_iodet", "arr_idet", "arr_samy", "arr_samz"]

        if method == "replace":
            index = tuple(
                [np.random.randint(0, shape[i], 10).tolist() for i in range(3)]
            )
            indexinfo = list(index)
            indexinfo.insert(self.stackdim, slice(None))
            info["index"] = tuple(indexinfo)
        elif method == "align":
            info["axes"] = (
                axis.factory(range(nstack)),
                axis.factory(units.Quantity(range(-nstack + 1, nhor), units="um")),
                axis.factory(units.Quantity(range(-nstack + 1, nvert), units="um")),
            )
        elif method == "expression":
            info["expression"] = "{}/{arr_iodet}"
            info["copy"] = ["arr_iodet"]
            info["select"] = nxprocess.results["counters"]["arr_iodet"]
        elif method == "copy":
            info["expression"] = "{}"
            info["copy"] = ["Fe-K", "Si-K"]
            info["skip"] = [s for s in signals if s not in info["copy"]] + counters
        elif method == "resample":
            info["encoders"] = {
                "y": {"counter": "arr_samy", "resolution": yencres},
                "z": {"counter": "arr_samz", "resolution": zencres},
            }
            info["shift"] = {"y": ypossmax, "z": zpossmax}

        groups = {}
        for group in range(2):
            groups["detector{:02d}".format(group)] = signals
        groups["counters"] = counters

        for group, signals in groups.items():
            group = nxprocess.results.nxdata(group).mkdir()
            _ = nxprocess.results.positioners()
            for name in signals:
                if name == "arr_samy":
                    data = yenc
                elif name == "arr_samz":
                    data = zenc
                else:
                    data = np.random.normal(size=shape)
                data = data.astype(dtype)
                if method == "crop":
                    data[:, 0, :] = np.nan
                    data[:, -2:, :] = np.nan
                    data[:, :, 0:2] = np.nan
                    data[:, :, -1] = np.nan
                    info["y"] = y[1][1:-2]
                    info["z"] = z[1][2:-1]
                elif method == "replace":
                    data[index] = -1
                elif method == "minlog":
                    mi = np.min(data) * 1.1
                    if mi == 0:
                        mi = -1
                    data -= mi
                elif method == "align":
                    hot = np.max(data) * 1.1
                    for i in range(nstack):
                        data[i, i, i] = hot
                group.add_signal(name, data=data)
            group.set_axes(x, y, z)

        try:
            yield nxprocess, info
        finally:
            # root.remove(recursive=True)
            pass

    def test_grid(self):
        with self._nxprocess() as proc:
            proc, info = proc
            grid = h5regulargrid.NXRegularGrid(proc)
            self._check_grid(grid)

            nxdata = proc.results["detector00"]
            grid = h5regulargrid.NXSignalRegularGrid(nxdata.signal)
            self._check_grid(grid)

    def test_crop(self):
        with self._nxprocess(method="crop") as proc1:
            proc1, info = proc1
            parameters = {
                "method": "crop",
                "default": "Si-K",
                "sliced": False,
                "reference": "Al-K",
                "nanval": np.nan,
            }
            proc2 = self._run_task(parameters, proc1)

            parameters["sliced"] = True
            proc3 = self._run_task(parameters, proc1)
            self._check_reproc(proc2, proc3)

            grid1 = h5regulargrid.NXRegularGrid(proc1)
            grid2 = h5regulargrid.NXRegularGrid(proc2)
            grid3 = h5regulargrid.NXRegularGrid(proc3)
            self.assertEqual(
                {sig.name for sig in grid1.signals}, {sig.name for sig in grid2.signals}
            )
            self.assertFalse(np.isnan(grid2.values).any())
            np.testing.assert_array_equal(grid2.values, grid3.values)

            for k, v in info.items():
                for ax in grid2.axes:
                    if ax.name == k:
                        np.testing.assert_array_equal(ax.magnitude, v)
                        break
                else:
                    assert False

            self.assertEqual(proc1.default.signal.name, parameters["default"])

    def test_replace(self):
        with self._nxprocess(method="replace") as proc1:
            proc1, info = proc1
            parameters = {
                "method": "replace",
                "default": "Si-K",
                "sliced": False,
                "old": -1,
                "new": -2,
            }
            proc2 = self._run_task(parameters, proc1)

            parameters["sliced"] = True
            proc3 = self._run_task(parameters, proc1)
            self._check_reproc(proc2, proc3)

            grid1 = h5regulargrid.NXRegularGrid(proc1)
            grid2 = h5regulargrid.NXRegularGrid(proc2)
            grid3 = h5regulargrid.NXRegularGrid(proc3)
            self.assertEqual(
                {sig.name for sig in grid1.signals}, {sig.name for sig in grid2.signals}
            )
            np.testing.assert_array_equal(grid2.values, grid3.values)
            np.testing.assert_array_equal(
                grid1.values[info["index"]], parameters["old"]
            )
            np.testing.assert_array_equal(
                grid2.values[info["index"]], parameters["new"]
            )

            self.assertEqual(proc1.default.signal.name, parameters["default"])

    def test_minlog(self):
        with self._nxprocess(method="minlog") as proc1:
            proc1, info = proc1
            parameters = {"method": "minlog", "sliced": False}
            proc2 = self._run_task(parameters, proc1)

            parameters["sliced"] = True
            proc3 = self._run_task(parameters, proc1)
            self._check_reproc(proc2, proc3)

            grid1 = h5regulargrid.NXRegularGrid(proc1)
            grid2 = h5regulargrid.NXRegularGrid(proc2)
            grid3 = h5regulargrid.NXRegularGrid(proc3)

            self.assertEqual(
                {sig.name for sig in grid1.signals}, {sig.name for sig in grid2.signals}
            )
            np.testing.assert_array_equal(grid2.values, grid3.values)
            np.testing.assert_array_equal(-np.log(grid1), grid3.values)

    def test_align(self):
        with self._nxprocess(method="align") as proc1:
            proc1, info = proc1
            parameters = {
                "method": "align",
                "alignmethod": "max",
                "reference": "Fe-K",
                "refimageindex": 0,
                "default": "Fe-K",
            }
            proc2 = self._run_task(parameters, proc1)

            grid2 = h5regulargrid.NXRegularGrid(proc2)
            axes = grid2.axes
            axes.pop(grid2.stackdim)
            for ax1, ax2 in zip(info["axes"], axes):
                self.assertEqual(ax1, ax2)

    def test_expression(self):
        with self._nxprocess(method="expression") as proc1:
            proc1, info = proc1
            copy = [{"method": "regex", "pattern": name} for name in info["copy"]]
            parameters = {
                "method": "expression",
                "expression": info["expression"],
                "copy": copy,
                "sliced": False,
            }
            proc2 = self._run_task(parameters, proc1)

            parameters["sliced"] = True
            proc3 = self._run_task(parameters, proc1)
            self._check_reproc(proc2, proc3)

            grid1 = h5regulargrid.NXRegularGrid(proc1)
            grid2 = h5regulargrid.NXRegularGrid(proc2)
            self._check_axes(grid1, grid2)

            index = grid1.locate(info["select"], None, None, None)
            norm = grid1[index]
            for i in range(grid1.shape[grid1.stackdim]):
                index = list(index)
                index[grid1.stackdim] = i
                index = tuple(index)
                data = grid1[index]
                if grid1.signals[i].name not in info["copy"]:
                    data = data / norm
                np.testing.assert_array_equal(data, grid2[index])

    def test_resample(self):
        with self._nxprocess(method="resample") as proc1:
            proc1, info = proc1

            params = ((["y"], ["y", "z"]), (True, False))

            for i, p in enumerate(itertools.product(*params), 1):
                axes, crop = p
                encoders = {k: v for k, v in info["encoders"].items() if k in axes}
                parameters = {
                    "name": "crop{}".format(i),
                    "method": "resample",
                    "encoders": encoders,
                    "crop": crop,
                }
                proc2 = self._run_task(parameters, proc1)

                grid1 = h5regulargrid.NXRegularGrid(proc1)
                grid2 = h5regulargrid.NXRegularGrid(proc2)

                # Check new axes position
                encoder_signals = {}
                offsets = proc2.results["encoder_offset"].read().tolist()
                offsets.insert(grid2.stackdim, 0)
                for ax1, ax2, offset in zip(grid1.axes, grid2.axes, offsets):
                    name = ax1.name
                    if name in axes:
                        n = int(np.ceil(info["shift"][name] / 2.0))
                        if crop:
                            ax1 = axis.factory(ax1[n:-n])
                        else:
                            add = np.arange(1, n + 1) * ax1.stepsize
                            addstart = ax1.start - add[::-1]
                            addend = ax1.end + add
                            x = (
                                addstart.magnitude.tolist()
                                + ax1.magnitude.tolist()
                                + addend.magnitude.tolist()
                            )
                            ax1 = axis.factory(units.Quantity(x, units=ax1.units))
                        resolution = units.Quantity(
                            parameters["encoders"][name]["resolution"],
                            units=1 / ax2.units,
                        )
                        encoder_signals[name] = (
                            (ax2.values * resolution + offset)
                            .to("dimensionless")
                            .magnitude
                        )
                    self._check_axis(ax1, ax2)

                # Check encoder signals
                signals = grid2.signal_names
                for axname, encinfo in encoders.items():
                    i = signals.index(encinfo["counter"])
                    enc = grid2[i, ...]
                    encnan = np.isnan(enc)
                    self.assertTrue(crop ^ encnan.any())

                    # Expected encoder values
                    encvalues = encoder_signals[axname]
                    index = [np.newaxis] * enc.ndim
                    if axname == "y":
                        index[1] = slice(None)
                    else:
                        index[2] = slice(None)
                    encvalues = encvalues[tuple(index)]

                    # Handle nan's
                    if not crop:
                        m = np.ones(enc.shape)
                        m[encnan] = np.nan
                        encvalues = encvalues * m
                        encvalues[encnan] = 999
                        enc[encnan] = 999

                    self.assertTrue(np.isclose(enc, encvalues).all())

    def test_copy(self):
        with self._nxprocess(method="copy") as proc1:
            proc1, info = proc1
            copy = [{"method": "regex", "pattern": name} for name in info["copy"]]
            skip = [{"method": "regex", "pattern": name} for name in info["skip"]]
            parameters = {
                "method": "expression",
                "expression": info["expression"],
                "copy": copy,
                "skip": skip,
                "sliced": False,
            }
            proc2 = self._run_task(parameters, proc1)

            parameters["sliced"] = True
            proc3 = self._run_task(parameters, proc1)
            self._check_reproc(proc2, proc3)

            grid1 = h5regulargrid.NXRegularGrid(proc2)
            grid2 = h5regulargrid.NXRegularGrid(proc2)
            signals1 = {s.name for s in grid1.signals if s.name not in info["skip"]}
            signals2 = {s.name for s in grid2.signals}
            self.assertEqual(signals1, signals2)

            index = [None] * grid2.ndim
            for s2 in grid2.signals:
                for s1 in grid1.signals:
                    if s1.name == s2.name and s1.parent.name == s2.parent.name:
                        break
                index[grid1.stackdim] = s1
                index1 = grid1.locate(*index)
                index[grid2.stackdim] = s2
                index2 = grid2.locate(*index)
                np.testing.assert_array_equal(grid1[index1], grid2[index2])

    @unittest.skip("h5py doesn't fully support concurrency")
    def test_concurrency(self):
        with self._nxprocess(method="copy") as proc1:
            proc1, info = proc1
            copy = [{"method": "regex", "pattern": name} for name in info["copy"]]
            skip = [{"method": "regex", "pattern": name} for name in info["skip"]]
            parameters = {
                "method": "expression",
                "expression": info["expression"],
                "copy": copy,
                "skip": skip,
                "sliced": False,
            }

            previoustask = utils.nxpathtotask(proc1)

            tasks = []
            threads = []
            nthreads = 5
            for i in range(nthreads):
                newtask = utils.create_task(dependencies=previoustask, **parameters)
                tasks.append(newtask)
                t = threading.Thread(target=newtask.run)
                threads.append(t)

            for t in threads:
                t.start()
            for t in threads:
                t.join()
            output = tasks[0].output
            for task in tasks:
                self.assertFalse(newtask.done)
                self.assertEqual(newtask.output, output)

            proc1.root.ls(recursive=True)

    def test_scenevis(self):
        with self._nxprocess(method="scenevis") as proc1:
            proc1, info = proc1
            outputparent = proc1.root["images"]
            obj1 = {
                "items": [
                    ("detector00/Fe-K", 0),
                    ("detector00/Si-K", 0),
                    ("detector00/Al-K", 0),
                ]
            }
            parameters = {
                "method": "scenevis",
                "outputparent": outputparent,
                "objects": [obj1],
                "instrument": "sxm",
                "title": "title1",
                "plot": False,
            }
        proc2 = self._run_task(parameters, proc1)
        parameters["title"] = "title2"
        proc3 = self._run_task(parameters, proc1)
        self._check_reproc(proc2, proc3)

    def _check_axes(self, grid1, grid2):
        for ax1, ax2 in zip(grid1.axes, grid2.axes):
            self._check_axis(ax1, ax2)

    def _check_axis(self, ax1, ax2):
        if ax1.type == "quantitative":
            self.assertEqual(ax1, ax2)
        else:
            self.assertEqual(len(ax1), len(ax2))
            for v1, v2 in zip(ax1, ax2):
                self.assertEqual(v1.name, v2.name)

    def _check_grid(self, grid):
        data = grid.values
        self.assertEqual(grid.shape, data.shape)
        self.assertEqual(grid.ndim, data.ndim)
        self.assertEqual(grid.size, data.size)
        np.testing.assert_array_equal(grid[:], data)

        indices = genindexing.genindexingn(
            data.shape, advanced=False, eco=False, nmax=50
        )
        for index in indices:
            np.testing.assert_array_equal(grid[index], data[index])

        for a, b in zip(grid, data):
            np.testing.assert_array_equal(a, b)


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_task_generic("test_grid"))
    testSuite.addTest(test_task_generic("test_copy"))
    testSuite.addTest(test_task_generic("test_concurrency"))
    testSuite.addTest(test_task_generic("test_crop"))
    testSuite.addTest(test_task_generic("test_replace"))
    testSuite.addTest(test_task_generic("test_minlog"))
    testSuite.addTest(test_task_generic("test_align"))
    testSuite.addTest(test_task_generic("test_expression"))
    testSuite.addTest(test_task_generic("test_resample"))
    testSuite.addTest(test_task_generic("test_scenevis"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

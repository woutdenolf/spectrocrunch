import unittest
import numpy as np
import os
import contextlib
import re
from testfixtures import TempDirectory

from .. import fluoxas
from ..run import run_sequential
from ...io import xiaedf
from ...process import nxresult
from ...utils import instance
from ...utils import listtools
from ...utils import cli
from ...materials import compoundfromname
from ...testutils.subtest import TestCase
from ...process.tests.xrfmap import XrfMapGenerator
from ...xrf.fit import OutputBuffer


logger = cli.getLogger(__name__, __file__)


class test_fluoxas(TestCase):
    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()

    @unittest.skipIf(compoundfromname.xraylib is None, "xraylib not installed")
    def test_process(self):
        self.xrfmap = XrfMapGenerator()
        self.xrfmap.generate(self.dir.path, "test")
        parameters = {
            "alignmethod": (None,),
            "cfgfileuse": (True, False),
            "include_detectors": self.xrfmap.include_detectors,
            "adddetectors": (True, False),
            "addbeforefit": (True, False),
            "quant": (True, False),
            "dtcor": (True, False),
            "stackdim": (0,),
            "correctspectra": (True, False),
        }
        _ = {
            "alignmethod": (None,),
            "cfgfileuse": (True,),
            "include_detectors": [(1, (0, 2))],
            "adddetectors": (False,),
            "addbeforefit": (False,),
            "quant": (True,),
            "dtcor": (True,),
            "stackdim": (0,),
            "correctspectra": (False,),
        }
        self.run_subtests(parameters, self._process)

    def _process(
        self,
        alignmethod=None,
        cfgfileuse=None,
        include_detectors=None,
        adddetectors=None,
        addbeforefit=None,
        quant=None,
        dtcor=None,
        stackdim=None,
        correctspectra=None,
    ):
        if not cfgfileuse and alignmethod is not None:
            return
        parameters = {}
        pymcaparams = {}
        prenormparams = {}
        postnormparams = {}
        alignparams = {}
        cropparams = {}
        replaceparams = {}
        commonparams = {}
        parameters["pymca"] = pymcaparams
        parameters["prealignnormalize"] = prenormparams
        parameters["postalignnormalize"] = postnormparams
        parameters["align"] = alignparams
        parameters["crop"] = cropparams
        parameters["replacenan"] = replaceparams
        parameters["common"] = commonparams
        commonparams["instrument"] = self.xrfmap.instrument
        commonparams["stackdim"] = stackdim

        nmaps, nlines, nspec, nchan, ndet = self.xrfmap.data.shape
        if include_detectors or include_detectors == 0:
            incdets_explicite = include_detectors
        else:
            incdets_explicite = tuple(range(ndet))

        alldetectors = tuple(listtools.flatten(incdets_explicite))
        adddetectorgroups = any(
            len(instance.asarray(dets)) > 1
            for dets in instance.asarray(incdets_explicite)
        )
        adddects_explicite = adddetectors and len(alldetectors) > 1
        addspectra = (adddects_explicite or adddetectorgroups) and addbeforefit
        fluxnormbefore = quant and correctspectra
        dtcorbefore = dtcor and (correctspectra or addspectra)
        newspectra = addspectra or fluxnormbefore or dtcorbefore
        if addspectra:
            if adddetectorgroups:
                seldetectors = [
                    tuple(instance.asarray(dets).tolist())
                    for dets in instance.asarray(incdets_explicite)
                ]
            else:
                seldetectors = [alldetectors]
        else:
            seldetectors = [(det,) for det in alldetectors]

        if cfgfileuse:
            cfgfiles = [self.xrfmap.xrfspectra[k]["cfgfile"] for k in seldetectors]
        else:
            cfgfiles = None

        if quant:
            geomparams = self.xrfmap.taskparams_geometry(seldetectors)
            parameters["geometry"] = geomparams
            pymcaparams.update(self.xrfmap.taskparams_pymca(seldetectors))
            prealignnormcounter = None
        else:
            geomparams = {}
            prealignnormcounter = "arr_norm"

        alignreference = None
        fitlabels = set(self.fitlabels(quant=quant))
        fitlabelsfile = fitlabels
        detcounterlabels = {"xmap_icr", "xmap_ocr", "xmap_x1c", "xmap_x2c"}
        counterlabels = {"arr_iodet", "arr_idet", "arr_norm"}
        calclabels = {
            "calc_transmission",
            "calc_absorbance",
            "calc_flux0",
            "calc_fluxt",
        }
        for label in fitlabels:
            if "Scatter" not in label:
                alignreference = label
                break
        refimageindex = 0
        usealign = alignmethod is not None and alignreference is not None

        # Data
        expectedgroups_data = []
        if addspectra:
            if adddetectorgroups:
                expectedgroups_data = {
                    "S{:d}".format(i + 1) for i in range(len(incdets_explicite))
                }
            else:
                if len(alldetectors) == 1:
                    expectedgroups_data = {"{:02d}".format(alldetectors[0])}
                else:
                    expectedgroups_data = {"S1"}
        else:
            expectedgroups_data = {
                "{:02d}".format(i) for i in list(listtools.flatten(incdets_explicite))
            }

        # Final groups
        if adddects_explicite:
            if adddetectorgroups:
                expectedgroups_result = ["S{:d}".format(len(incdets_explicite) + 1)]
            else:
                expectedgroups_result = ["S1"]
        elif adddetectorgroups:
            expectedgroups_result = [
                "S{:d}".format(i + 1) for i in range(len(incdets_explicite))
            ]
        else:
            expectedgroups_result = ["{:02d}".format(i) for i in alldetectors]

        if alignreference:
            alignreference = "/detector{}/{}".format(
                expectedgroups_result[0], alignreference
            )
        expectedgroups_result = ["counters"] + [
            "detector" + det for det in expectedgroups_result
        ]
        expectedgroups_result = set(expectedgroups_result)

        # Processes
        expected_nxprocess = ["pymca.1"]
        if prealignnormcounter is not None:
            expected_nxprocess.append("normalize.1")
        if usealign:
            expected_nxprocess.append("align.1")
            expected_nxprocess.append("crop.1")

        with self._destpath_context() as destpath:
            radix = self.xrfmap.radix
            commonparams["outputparent"] = os.path.join(
                destpath.path, radix + ".h5::/" + radix
            )
            geomparams["outputparent"] = os.path.join(
                destpath.path, radix + ".h5::/" + radix
            )
            pymcaparams["sourcepaths"] = [self.xrfmap.path]
            pymcaparams["scannames"] = [radix]
            pymcaparams["scannumbers"] = [self.xrfmap.scannumbers]
            pymcaparams["pymcacfg"] = cfgfiles
            pymcaparams["dtcor"] = dtcor
            pymcaparams["adddetectors"] = adddetectors
            pymcaparams["addbeforefit"] = addbeforefit
            pymcaparams["correctspectra"] = correctspectra
            pymcaparams["include_detectors"] = include_detectors
            pymcaparams["counters"] = ["arr_norm"]
            pymcaparams["fastfitting"] = True
            prenormparams["counter"] = prealignnormcounter
            alignparams["alignmethod"] = alignmethod
            alignparams["reference"] = alignreference
            alignparams["refimageindex"] = refimageindex
            alignparams["plot"] = False
            cropparams["crop"] = usealign
            replaceparams["replacenan"] = False

            for repeat in range(2):
                tasks = fluoxas.tasks(**parameters)
                if repeat:
                    for task in tasks:
                        self.assertTrue(task.done)
                    continue
                else:
                    for task in tasks:
                        self.assertFalse(task.done)
                    run_sequential(tasks)
                    for task in tasks:
                        self.assertTrue(task.done)
                    nxprocess = tasks[-1].output

                # Check generated spectra (files)
                if newspectra:
                    corlabel = ""
                    if dtcorbefore:
                        corlabel = corlabel + "dt"
                    if fluxnormbefore:
                        corlabel = corlabel + "fl"
                    if corlabel:
                        radixout = "{}_{}cor".format(radix, corlabel)
                    else:
                        radixout = radix

                    if addspectra:
                        expected = [
                            "{}_xia{}_{:04d}_0000_{:04d}.edf".format(
                                radixout, det, mapnum, linenum
                            )
                            for det in expectedgroups_data
                            for mapnum in range(nmaps)
                            for linenum in range(nlines)
                        ]
                    else:
                        expected = [
                            "{}_xia{}_{:04d}_0000_{:04d}.edf".format(
                                radixout, det, mapnum, linenum
                            )
                            for det in expectedgroups_data
                            for mapnum in range(nmaps)
                            for linenum in range(nlines)
                        ]
                    xrfspectra_subdir = os.path.join(
                        "{}_pymca.1".format(radix), "xrfspectra"
                    )
                    destpath.compare(
                        sorted(expected),
                        path=xrfspectra_subdir,
                        files_only=True,
                        recursive=False,
                    )
                else:
                    radixout = radix

                # Check pymca fit output (files)
                if cfgfileuse:
                    if addspectra:
                        expected = [
                            "{}_xia{}_{:04d}_0000_{}.edf".format(
                                radixout, det, mapnum, label
                            )
                            for det in expectedgroups_data
                            for mapnum in range(nmaps)
                            for label in fitlabelsfile
                        ]
                        expected.extend(
                            [
                                "{}_xia{}_{:04d}_0000.cfg".format(radixout, det, mapnum)
                                for det in expectedgroups_data
                                for mapnum in range(nmaps)
                            ]
                        )
                    else:
                        expected = [
                            "{}_xia{}_{:04d}_0000_{}.edf".format(
                                radixout, det, mapnum, label
                            )
                            for det in expectedgroups_data
                            for mapnum in range(nmaps)
                            for label in fitlabelsfile
                        ]
                        expected.extend(
                            [
                                "{}_xia{}_{:04d}_0000.cfg".format(radixout, det, mapnum)
                                for det in expectedgroups_data
                                for mapnum in range(nmaps)
                            ]
                        )
                    fitresults_subdir = os.path.join(
                        "{}_pymca.1".format(radix), "pymcaresults"
                    )
                    if OutputBuffer is None:
                        destpath.compare(
                            sorted(expected),
                            path=fitresults_subdir,
                            files_only=True,
                            recursive=False,
                        )
                    else:
                        expected = set(expected)
                        actual = set(
                            os.listdir(os.path.join(destpath.path, fitresults_subdir))
                        )
                        self.assertEqual(actual & expected, expected)

                # Check top-level output directory (h5 files)
                expected = []
                if cfgfiles or newspectra:
                    expected.append("{}_pymca.1".format(radix))
                h5file = "{}.h5".format(radix)
                expected.append(h5file)
                destpath.compare(sorted(expected), files_only=True, recursive=False)

                # Check NXprocess groups
                entry = nxprocess.nxentry()
                for name in expected_nxprocess:
                    self.assertTrue(name in entry)
                self.assertEqual(nxprocess.name, expected_nxprocess[-1])

                # Check NXdata groups
                groups, axes, stackdim = nxresult.regulargriddata(nxprocess)
                self.assertEqual(set(groups.keys()), expectedgroups_result)
                for group, signals in groups.items():
                    if group == "counters":
                        if quant:
                            expectedsubgroups = counterlabels | calclabels
                        else:
                            expectedsubgroups = counterlabels
                    else:
                        if cfgfileuse:
                            expectedsubgroups = detcounterlabels | fitlabels
                        else:
                            expectedsubgroups = detcounterlabels
                    self.assertEqual({sig.name for sig in signals}, expectedsubgroups)

                # Check generated spectra (data)
                if newspectra:
                    # Apply DT correction
                    if dtcorbefore:
                        data0 = (
                            self.xrfmap.stats[..., xiaedf.xiadata.STICR, :]
                            / self.xrfmap.stats[..., xiaedf.xiadata.STOCR, :]
                        )
                        data0 = self.xrfmap.data * data0[..., np.newaxis, :]
                    else:
                        data0 = self.xrfmap.data.copy()

                    # Apply flux normalization
                    if fluxnormbefore:
                        data0 /= self.xrfmap.ctrs["arr_norm"][
                            ..., np.newaxis, np.newaxis
                        ]

                    # Add spectra
                    if addspectra:
                        if adddetectorgroups:
                            data0 = np.stack(
                                [
                                    data0[..., instance.asarray(ind)].sum(axis=-1)
                                    for ind in incdets_explicite
                                ],
                                axis=-1,
                            )
                        else:
                            data0 = data0[..., alldetectors].sum(axis=-1)[
                                ..., np.newaxis
                            ]
                    else:
                        data0 = data0[..., tuple(sorted(alldetectors))]

                    # Saved spectra
                    stack = xiaedf.xiastack_radix(
                        os.path.join(destpath.path, xrfspectra_subdir), radixout
                    )
                    data2 = stack.data

                    # Check spectra are equal
                    np.testing.assert_allclose(data0, data2, rtol=1e-6)

                # Check fit results
                if cfgfileuse and dtcor:
                    for group, signals in groups.items():
                        if not group.isdetector:
                            continue
                        if group.issum:
                            if adddetectorgroups:
                                if group.number > len(incdets_explicite):
                                    dets = tuple(sorted(alldetectors))
                                else:
                                    dets = tuple(
                                        instance.asarray(
                                            incdets_explicite[group.number - 1]
                                        ).tolist()
                                    )
                            else:
                                dets = tuple(sorted(alldetectors))
                        else:
                            dets = (group.number,)
                        logger.debug(
                            "Check fit result for sum of xrfdata {}".format(dets)
                        )
                        info = self.xrfmap.xrfspectra[dets]
                        for signal in signals:
                            if "xmap" in signal.name:
                                continue
                            dataset = signal.read()
                            if stackdim == 1:
                                grpdata = np.moveaxis(dataset, 1, 0)
                            elif stackdim == 2:
                                grpdata = np.moveaxis(dataset, 2, 0)
                            else:
                                grpdata = dataset[:]
                            self._assert_fitresult(signal.name, grpdata, info)
            # repeat
        # destpath context

    def fitlabels(self, quant=False):
        labels = []
        if OutputBuffer is not None:
            labels += ["Constant"]
        labels += list(self.xrfmap.labels)
        if quant:
            labels += [
                "w" + label
                for label in labels
                if "Scatter" not in label and label != "Constant"
            ]
        labels = [label.replace("-", "_") for label in labels]
        return labels

    @contextlib.contextmanager
    def _destpath_context(self):
        destpath = TempDirectory()
        yield destpath
        destpath.cleanup()

    def _assert_fitresult(self, grpname, grpdata, info):
        if "Scatter" in grpname or "chisq" in grpname or "Constant" in grpname:
            return
        grpname = str(grpname).replace("_", "-")
        m = re.match("Scatter-(Compton|Peak)([0-9]+)", grpname)
        if m:
            grpname = m.group(1)
            if grpname == "Peak":
                grpname = "Rayleigh"
            values1 = [
                peakareas[0][grpname][int(m.group(2))]
                for peakareas in info["peakareas"]
            ]
            values2 = [
                peakareas[1][grpname][int(m.group(2))]
                for peakareas in info["peakareas"]
            ]
        else:
            if grpname.startswith("w"):
                grpname = grpname[1:]
                values1 = [
                    massfractions[0][grpname] for massfractions in info["massfractions"]
                ]
                values2 = [
                    massfractions[1][grpname] for massfractions in info["massfractions"]
                ]
            else:
                values1 = [peakareas[0][grpname] for peakareas in info["peakareas"]]
                values2 = [peakareas[1][grpname] for peakareas in info["peakareas"]]
        for data, v1, v2 in zip(grpdata, values1, values2):
            mask = data == np.nanmax(data)
            np.testing.assert_allclose(data[~mask], v1, rtol=1e-4, err_msg=grpname)
            np.testing.assert_allclose(data[mask], v2, rtol=1e-4, err_msg=grpname)

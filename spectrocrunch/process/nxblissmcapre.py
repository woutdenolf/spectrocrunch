# -*- coding: utf-8 -*-

import re
from . import nxprocess
from ..io import nxfs
from ..math.slicing import slice_generator
import numpy


class Task(nxprocess.Task):
    """Prepare raw data for XRF fitting"""

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()

        self.required_parameters |= {
            "detectors",  # NXdetector names of MCA's
            "preset_time",
        }
        self.optional_parameters |= {"dtcor", "add", "detector_suffix", "counters"}

        parameters = self.parameters
        parameters["dtcor"] = parameters.get("dtcor", True)
        parameters["add"] = parameters.get("add", True)
        parameters["detector_suffix"] = parameters.get("detector_suffix", "cor")
        parameters["counters"] = parameters.get("counters", [])

    def _execute(self):
        if len(self.dependencies) != 1:
            raise RuntimeError("Expected 1 dependency")
        source_nxentry = self.dependencies[0]
        if source_nxentry.is_not_nxclass(u"NXentry"):
            raise RuntimeError("Dependency needs to be an NXentry")
        dest_nxentry = self.outputparent

        if source_nxentry.device == dest_nxentry.device:
            mode = "a"
        else:
            mode = "r"

        parameters = self.parameters
        with source_nxentry.open(mode=mode):
            detectors = list(
                self._matching_names(
                    source_nxentry["instrument"], self.parameters["detectors"]
                )
            )
            if not detectors:
                raise RuntimeError("No matching detector names found")

            with dest_nxentry.open():
                # Add NXdetector's
                if parameters["add"]:
                    assert parameters[
                        "dtcor"
                    ], "Cannot add MCA's without deadtime correction"
                    self._mca_add(source_nxentry["instrument"], detectors, "data")
                elif parameters["dtcor"]:
                    self._mca_dtcor(source_nxentry["instrument"], detectors)
                self._add_counters(source_nxentry["instrument"])

    def _matching_names(self, path, patterns):
        if not patterns:
            return
        names = [f.name for f in path.listdir()]
        for pattern in patterns:
            pattern = re.compile(pattern)
            for name in names:
                if pattern.match(name):
                    yield name

    def _mca_add(self, nxinstrument, detectors, dsetname, outname=None):
        datasets = []
        livetimes = []
        preset_time = self.parameters["preset_time"]
        for source_name in detectors:
            det = nxinstrument[source_name]
            with det[dsetname].open() as dset:
                shape = dset.shape
                dtype = dset.dtype
                datasets.append(dset)

            detlt = None
            for name in ["live_time", "output_live_time"]:
                detlt = det[name]
                if detlt.exists:
                    break

            with detlt.open() as dset:
                livetimes.append(dset)
                dtype = (numpy.array(0, dset.dtype) * preset_time).dtype

        mca = self._mca_nxdetector("mca" + self.parameters["detector_suffix"])
        it = slice_generator(shape[::-1], dtype)
        if not outname:
            outname = dsetname
        with mca[outname].open(shape=shape, dtype=dtype) as result:
            for idx in it:
                idx = idx[::-1]
                idxlt = idx[:-1]
                result[idx] = (
                    datasets[0][idx][()]
                    * preset_time
                    / livetimes[0][idxlt][..., numpy.newaxis]
                )
                for dset, lt in zip(datasets[1:], livetimes[1:]):
                    result[idx] += (
                        dset[idx][()] * preset_time / lt[idxlt][..., numpy.newaxis]
                    )
            mca["preset_time"].write(data=preset_time)

    def _mca_dtcor(self, nxinstrument, detectors):
        preset_time = self.parameters["preset_time"]
        for source_name in detectors:
            det = nxinstrument[source_name]
            with det["data"].open() as dset:

                detlt = None
                for name in ["live_time", "output_live_time"]:
                    detlt = det[name]
                    if detlt.exists:
                        break

                with detlt.open() as lt:
                    shape = dset.shape
                    dtype = (numpy.array(0, dset.dtype) * preset_time).dtype
                    it = slice_generator(shape[::-1], dtype)
                    mca = self._mca_nxdetector(
                        source_name + self.parameters["detector_suffix"]
                    )
                    with mca["data"].open(shape=shape, dtype=dtype) as result:
                        for idx in it:
                            idx = idx[::-1]
                            idxlt = idx[:-1]
                            result[idx] = (
                                dset[idx] * preset_time / lt[idxlt][..., numpy.newaxis]
                            )
                    mca["preset_time"].write(data=preset_time)

    def _mca_nxdetector(self, name):
        mca = self.temp_nxresults.nxdetector(name)
        mca["type"].write(data="mca")
        self.temp_nxresults[name].link(mca)
        return mca

    def _add_counters(self, nxinstrument):
        dest = self.temp_nxresults.nxinstrument()
        for name in self.parameters["counters"]:
            dest[name].link(nxinstrument[name])
            # Softlinks cause problems later on (softlink to external links)
            nxinstrument[name].copy(
                self.temp_nxresults[name], follow=True, dereference=True
            )
            # self.temp_nxresults[name].link(dest[name])

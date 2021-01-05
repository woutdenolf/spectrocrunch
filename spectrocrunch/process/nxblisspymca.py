# -*- coding: utf-8 -*-

from . import nxprocess
from ..utils import instance
from ..io import fs
from ..io import nxfs
from ..xrf.fit import PerformBatchFitHDF5

from PyMca5.PyMcaIO import ConfigDict


class Task(nxprocess.Task):
    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()

        self.required_parameters |= {"pymcacfg"}
        self.optional_parameters |= {
            "fluxcounter",
            "transmissioncounter",
            "counters",
            "fastfitting",
            "mlines",
            "addhigh",
            "diagnostics",
        }

        parameters = self.parameters

        pymcacfg = parameters.get("pymcacfg", None)
        if pymcacfg is None:
            pymcacfg = []
        elif isinstance(pymcacfg, str):
            pymcacfg = [pymcacfg]
        parameters["pymcacfg"] = [
            ConfigDict.ConfigDict(filelist=cfg) if instance.isstring(cfg) else cfg
            for cfg in pymcacfg
        ]

        parameters["fluxcounter"] = parameters.get("fluxcounter", None)
        parameters["transmissioncounter"] = parameters.get("transmissioncounter", None)
        parameters["fastfitting"] = parameters.get("fastfitting", True)
        parameters["mlines"] = parameters.get("mlines", {})
        parameters["addhigh"] = parameters.get("addhigh", 0)
        parameters["diagnostics"] = parameters.get("diagnostics", False)

    def _execute(self):
        self._outuris = []
        self._detnames = []
        if len(self.dependencies) != 1:
            raise RuntimeError("Expected 1 blissmcapre task")
        if self.dependencies[0].method != "blissmcapre":
            raise RuntimeError("Expected 1 blissmcapre task")
        parameters = self.parameters

        cfgs = parameters["pymcacfg"]
        detectors = list(self.detectors)
        if len(detectors) != len(cfgs):
            if len(cfgs) != 1:
                raise RuntimeError(
                    "Provide 1 or {} pymca config files".format(len(detectors))
                )
            cfgs = cfgs * len(detectors)

        # Pymca fitting
        nxentry = self.temp_nxprocess.parent
        for cfg, detector in zip(cfgs, detectors):
            outuri = nxentry[self.temp_outputname + ":" + detector.name]
            self._outuris.append(outuri)
            self._detnames.append(detector.name)
            PerformBatchFitHDF5(
                [str(detector["data"])],
                cfg,
                outuri,
                energy=None,
                mlines=parameters["mlines"],
                quant=None,
                fast=parameters["fastfitting"],
                addhigh=parameters["addhigh"],
                diagnostics=parameters["diagnostics"],
            )

        nxentry["plotselect"].remove()

        # Add counters to pymca results
        counters = list(self.counters)
        if counters:
            for outuri in self._outuris:
                dest = outuri["results"].nxdata("counters")
                for counter in counters:
                    dest.add_signal(name=counter.name, path=counter["data"])

        # Add pymca results
        nxdata_names = ["parameters", "concentrations"]
        if counters:
            nxdata_names.append("counters")
        with nxentry.open():
            results = self.temp_nxresults
            outuri = self._outuris[0]
            # Copy results of one detector
            for nxdataname in nxdata_names:
                nxdatain = outuri["results"][nxdataname]
                if not nxdatain.exists:
                    continue
                nxdataout = results[nxdataname]
                nxdatain.copy(nxdataout, follow=True, dereference=True)
                # for path in nxdataout:
                #    if path.islink:
                #        path.remove()
            # Add results of other detectors
            for outuri in self._outuris[1:]:
                for nxdataname in nxdata_names:
                    nxdatain = outuri["results"][nxdataname]
                    if not nxdatain.exists:
                        continue
                    nxdataout = results[nxdataname]
                    for singalin in nxdatain.signals:
                        signalout = nxdataout[singalin.name]
                        try:
                            with signalout.open() as dsetout:
                                with singalin.open() as dsetin:
                                    dsetout[()] += dsetin[()]
                        except Exception:
                            pass

    def _remove_temp_output(self):
        super(Task, self)._remove_temp_output()
        if not hasattr(self, "_outuris"):
            return
        for outuri in self._outuris:
            outuri.remove(recursive=True)

    def _rename_temp_output(self):
        super(Task, self)._rename_temp_output()
        if not hasattr(self, "_outuris"):
            return
        outparent = self.outputparent
        outname = self.outputname
        for old, detname in zip(self._outuris, self._detnames):
            new = outparent[outname + "." + detname]
            try:
                old.rename(new)
            except fs.AlreadyExists:
                if new.exists:
                    old.remove(recursive=True)

    @property
    def counters(self):
        for detector in self.previous_outputs[0].results:
            try:
                if detector["type"].read() != "mca":
                    yield detector
            except fs.Missing:
                yield detector

    @property
    def detectors(self):
        for detector in self.previous_outputs[0].results:
            try:
                if detector["type"].read() == "mca":
                    yield detector
            except fs.Missing:
                pass

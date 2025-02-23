from . import nxqxrf_dependent
from . import nxprocess
from . import nxresult
from . import nxlazy
from . import axis
from ..io import spec
from ..io import xiaedf
from ..utils import listtools
from ..utils import instance
from ..utils import units
from ..xrf.fit import PerformBatchFit

import numpy as np
from collections import OrderedDict
import logging
from PyMca5.PyMcaIO import ConfigDict


logger = logging.getLogger(__name__)


class Task(nxqxrf_dependent.Task, nxprocess.Task):

    DEFAULT_STACKDIM = 0

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()

        self.required_parameters |= {
            # Input
            "sourcepaths",
            "counter_reldir",
            "scannames",
            "scannumbers",
            "counters",
            "fluxcounter",
            "transmissioncounter",
            # Meta data
            "metadata",
            "edfheader",
            "units",
            # Data correction
            "dtcor",
            "correctspectra",
            "adddetectors",
            "addbeforefit",
            # Configuration for fitting
            "pymcacfg",
            "mlines",
            "addhigh",
            "fastfitting",
            "exclude_detectors",
            "include_detectors",
        }

        self.optional_parameters |= {
            "quantification",
            "stackdim",
            "ignore_energy",  # Do not use the energy for fitting (use the one from pymca)
        }

        parameters = self.parameters
        parameters["stackdim"] = parameters.get("stackdim", self.DEFAULT_STACKDIM)
        parameters["quantification"] = parameters.get("quantification", None)
        parameters["ignore_energy"] = parameters.get("ignore_energy", False)

        edfheader = parameters["edfheader"]
        edffields = (
            "speclabel",
            "slowlabel",
            "fastlabel",
            "stackvalue",
            "time",
            "axesnamemap",
            "compensationmotors",
        )
        defaults = (None, None, None, None, None, {}, {})

        parameters["edfheader"] = {
            k: edfheader.get(k, default) for k, default in zip(edffields, defaults)
        }

        pymcacfg = parameters.get("pymcacfg", None)
        if pymcacfg is None:
            pymcacfg = []
        else:
            pymcacfg = [
                ConfigDict.ConfigDict(filelist=cfg) if instance.isstring(cfg) else cfg
                for cfg in pymcacfg
            ]
        parameters["pymcacfg"] = pymcacfg

    def _execute(self):
        self._preparestacks()
        self._processstacks()

    def _preparestacks(self):
        self._prepare_input()
        self._prepare_output()
        self._prepare_adddetector()
        self._prepare_fluxnorm()
        self._prepare_dtcor()
        self._stack_add_counters()
        self._prepare_xrfquant()
        self._process_xiastackraw()
        self._stack_add_flux()
        self._stack_add_xrffit()
        self._postcorrect()
        self._sort_stacks()

    def _processstacks(self):
        # NXdata positioners
        positioners = self.temp_nxresults.positioners()
        for ax in self.axes.values():
            positioners.add_axis(ax.name, ax.values, title=ax.title)

        # Processing info
        info = self.temp_nxresults.nxcollection("info")
        for k, v in self.procinfo.items():
            info[k].mkfile(data=v)

        # Processing axes
        positioners = self.temp_nxresults.nxcollection("stackaxes")
        for ax in self.infoaxes.values():
            positioners.add_axis(ax.name, ax.values, title=ax.title)

        # Processing datasets
        self._exportgroups()

    def _prepare_input(self):
        # Check data
        npaths = len(self.parameters["sourcepaths"])
        if npaths != len(self.parameters["scannames"]):
            raise ValueError(
                "Number of scan names must be the same as number of source paths."
            )
        if npaths != len(self.parameters["scannumbers"]):
            raise ValueError(
                "Number of scan numbers must be the same as number of source paths."
            )

        # Get image stack
        self.xiastackraw = xiaedf.xiastack_mapnumbers(
            self.parameters["sourcepaths"],
            self.parameters["scannames"],
            self.parameters["scannumbers"],
        )
        if self.xiastackraw.isempty:
            raise IOError(
                "Cannot find data: {}".format(self.xiastackraw.filedescription)
            )

        dshape = self.xiastackraw.dshape  # nstack, nrow, ncol, nchan, ndet
        self.ndetorg = dshape[-1]  # xia00, xia01, xiaS0, ...
        self.nstack = dshape[0]

    def _prepare_output(self):
        self.stacks = OrderedDict()
        self.axes_names = [None] * 3
        self.axes = {}
        self.procinfo = {}
        self.infoaxes = {}

        self.outstackdim = self.parameters["stackdim"]
        if self.outstackdim == 0:
            self.outimgdim = [1, 2]
        elif self.outstackdim == 1:
            self.outimgdim = [0, 2]
        else:
            self.outimgdim = [0, 1]

    @property
    def outdatapath(self):
        return self.temp_localpath["xrfspectra"]

    @property
    def outfitpath(self):
        return self.temp_localpath["pymcaresults"]

    def _prepare_adddetector(self):
        # include_detectors = [1, (0, 2)]
        # When adding detector spectra, this means we will fit two
        # detector groups S1 and S2, with potentially two config files
        include_detectors = self.parameters["include_detectors"]
        used_detectors = list(listtools.flatten(include_detectors))

        # Detector include/exclude
        self.xiastackraw.exclude_detectors = self.parameters["exclude_detectors"]
        self.xiastackraw.include_detectors = used_detectors
        tmp = self.xiastackraw.detectors_used  # e.g. ['00', '01', 'S0']
        if used_detectors:
            # Keep order provided by 'include_detectors'
            tmp = [int(det) if det.isdigit() else det for det in tmp]
            used_detectors = [det for det in used_detectors if det in tmp]
        else:
            used_detectors = list(sorted(tmp))

        # Detector groups for fitting
        if any(
            len(listtools.aslist(lst)) > 1
            for lst in listtools.aslist(include_detectors)
        ):
            self._detector_groups = OrderedDict(
                (nxresult.Group("S{:d}".format(i)), listtools.aslist(singledets))
                for i, singledets in enumerate(include_detectors, 1)
            )
        else:
            tmp = [nxresult.Group(det) for det in used_detectors]
            tmp = [det.number for det in tmp if det.issum]
            if tmp:
                num = max(tmp) + 1
            else:
                num = 1
            self._detector_groups = {
                nxresult.Group("S{:d}".format(num)): used_detectors
            }

        # Do we need to add all detectors?
        adddetectors = (len(used_detectors) > 1) and self.parameters["adddetectors"]

        # Do we need to add detector groups?
        adddetectorgroups = len(self._detector_groups) > 1

        # Do we need to add spectra (all or in groups)?
        self.addspectra = (adddetectors or adddetectorgroups) and self.parameters[
            "addbeforefit"
        ]
        self.xiastackraw.detectorsum(self.addspectra)

        # List of detectors to fit
        if self.addspectra:
            self._detectors_to_fit = list(self._detector_groups.keys())
        else:
            self._detectors_to_fit = [nxresult.Group(det) for det in used_detectors]
        self.ndetfit = len(self._detectors_to_fit)
        logger.info(
            "Fitted detectors (in order): {}".format(
                [det.xialabel for det in self._detectors_to_fit]
            )
        )

        # Sum after fitting:
        self._detectors_sumto = OrderedDict()
        if adddetectorgroups or adddetectors:
            # Add single detectors to their corresponding subgroup
            for dest, singledets in self._detector_groups.items():
                for det in singledets:
                    self._detectors_sumto[nxresult.Group(det)] = dest
        if adddetectorgroups and adddetectors:
            # Add subgroups to the sumgroup
            num = max(det.number for det in self._detectors_sumto)
            dest = nxresult.Group("S{:d}".format(num + 1))
            for det in self._detector_groups:
                self._detectors_sumto[nxresult.Group(det)] = dest

    def _prepare_dtcor(self):
        dtcor = self.parameters["dtcor"]
        self.dtcorbefore = dtcor and (
            self.parameters["correctspectra"] or self.addspectra
        )
        self.dtcorafter = dtcor and not self.dtcorbefore
        self.procinfo["dtneeded"] = False
        self.xiastackraw.dtcor(self.dtcorbefore)

    def _prepare_fluxnorm(self):
        self.fluxnorm = self.qxrfgeometry is not None
        self.fluxnormbefore = self.fluxnorm and self.parameters["correctspectra"]
        self.fluxnormafter = self.fluxnorm and not self.fluxnormbefore
        self.procinfo["fluxnorm"] = self.fluxnorm

    @property
    def units(self):
        return self.parameters["units"]

    def _prepare_xrfquant(self):
        if not self.fluxnorm:
            return
        qxrfgeometry = self.qxrfgeometry
        ngeometries = len(qxrfgeometry.xrfgeometries)
        if ngeometries != self.ndetfit:
            raise RuntimeError(
                "You need {} detector geometries, {} provides.".format(
                    self.ndetfit, ngeometries
                )
            )

        def tile(x):
            return np.tile(np.array(x)[np.newaxis, :], (self.nstack, 1))

        self._add_info_axis("refflux", defaultunits="Hz")
        self._add_info_axis("refexpotime", defaultunits="s")
        self._add_info_axis(
            "activearea",
            values=tile(qxrfgeometry.xrf_activeareas.to("cm**2").magnitude),
            defaultunits="cm**2",
        )
        self._add_info_axis(
            "anglein", values=tile(qxrfgeometry.xrf_anglesin), defaultunits="deg"
        )
        self._add_info_axis(
            "angleout", values=tile(qxrfgeometry.xrf_anglesout), defaultunits="deg"
        )
        self._add_info_axis(
            "sampledetdistance",
            values=np.full((self.nstack, self.ndetfit), np.nan, dtype=np.float32),
            defaultunits="cm",
        )

        if self.fluxnormbefore:
            self._add_info_axis("i0_to_norm_offset")
            self._add_info_axis("i0_to_norm_factor")

        for imageindex, xiaimage in enumerate(self.xiastackraw):
            energy = (
                self.axes[self.axes_names[self.outstackdim]][imageindex]
                .to("keV")
                .magnitude
            )
            if not np.isnan(energy):
                time = self.infoaxes["expotime"][imageindex]
                if np.isnan(time.magnitude):
                    time = None

                (
                    xrfnormop,
                    self.infoaxes["refflux"][imageindex],
                    self.infoaxes["refexpotime"][imageindex],
                    self.infoaxes["expotime"][imageindex],
                ) = qxrfgeometry.xrfnormop(energy, expotime=time)

                if self.fluxnormbefore:
                    xiaimage.localnorm(self.parameters["fluxcounter"], func=xrfnormop)
                    self.infoaxes["i0_to_norm_offset"][imageindex] = xrfnormop.b
                    self.infoaxes["i0_to_norm_factor"][imageindex] = xrfnormop.m

                pos = self.infoaxes["xrfdetectorposition"][imageindex, :]
                if np.isfinite(pos).all():
                    qxrfgeometry.xrf_positions = pos
                self.infoaxes["xrfdetectorposition"][imageindex, :] = (
                    qxrfgeometry.xrf_positions.to("cm").magnitude
                )
                self.infoaxes["sampledetdistance"][imageindex, :] = (
                    qxrfgeometry.xrf_distances.to("cm").magnitude
                )

    def _process_xiastackraw(self):
        self.xiastackraw.exclude_detectors = self.parameters["exclude_detectors"]
        self.xiastackraw.include_detectors = list(
            listtools.flatten(self.parameters["include_detectors"])
        )
        msg = "Corrections before XRF fitting:"
        if self.fluxnormbefore:
            msg += " flux normalization"
        if self.dtcorbefore:
            if not msg.endswith(":"):
                msg += ","
            msg += " deadtime correction"
        if self.addspectra:
            if not msg.endswith(":"):
                msg += ","
            msg += " add detectors"
        if msg.endswith(":"):
            msg += " none (fit raw XRF spectra)"
        logger.info(msg)

        if self.dtcorbefore or self.addspectra or self.fluxnormbefore:
            label = ""
            if self.dtcorbefore:
                label = label + "dt"
            if self.fluxnormbefore:
                label = label + "fl"
            if label:
                label = label + "cor"
                radix = [
                    "{}_{}".format(radix, label)
                    for radix in self.parameters["scannames"]
                ]
            else:
                radix = self.parameters["scannames"]

            # Raises error when it already exists
            outdatapath = self.outdatapath
            outdatapath.mkdir(force=False)

            self.xiastackproc = xiaedf.xiastack_mapnumbers(
                outdatapath.path, radix, self.parameters["scannumbers"]
            )

            # Handle existing data: will not be the case as long as we keep outdatapath.mkdir(force=False)
            shape = self.xiastackproc.dshape
            if shape:
                create = shape[:-1] != self.xiastackraw.dshape[:-1]
            else:
                create = True
            if create:
                # Only spectra, not counters (more efficient that way)
                self.xiastackproc.overwrite(True)
                if self.addspectra:
                    logger.info(
                        "Creating corrected XRF spectra: {}".format(
                            list(self._detector_groups.keys())
                        )
                    )
                    for det, singledets in self._detector_groups.items():
                        self.xiastackraw.include_detectors = singledets
                        self.xiastackproc.save(
                            self.xiastackraw, xialabels=[det.xialabel]
                        )
                else:
                    xialabels = self.xiastackraw.xialabels_used
                    logger.info("Creating corrected XRF spectra: {}".format(xialabels))
                    self.xiastackproc.save(self.xiastackraw, xialabels=xialabels)
            else:
                logger.info("Corrected XRF spectra already exist")
        else:
            self.xiastackproc = self.xiastackraw

    def _stack_add_counters(self):
        self.xiastackraw.exclude_detectors = self.parameters["exclude_detectors"]
        self.xiastackraw.include_detectors = list(
            listtools.flatten(self.parameters["include_detectors"])
        )

        # Counter directory relative to the XIA files
        self.xiastackraw.counter_reldir(self.parameters["counter_reldir"])

        # Check counters
        countersfound = set(self.xiastackraw.counterbasenames())
        counters = countersfound.intersection(self.parameters["counters"])
        if self.parameters["metadata"] == "xia":
            metacounters = "xia"
        else:
            if countersfound:
                metacounters = next(iter(countersfound))
            else:
                logger.warning(
                    "Metacounters for {} are not found".format(self.xiastackraw)
                )
                metacounters = None

        logger.info("Processing counters: {}".format(counters))

        # Extract metadata and counters from raw stack
        self.counters = set()
        for imageindex, xiaimage in enumerate(self.xiastackraw):
            h = xiaimage.header(source=metacounters)
            parsedheader = self._getscanparameters(h)

            # Prepare axes and stackinfo
            if imageindex == 0:
                # Image axes: may be different for each image due to drift compenation,
                #             but keep the values of the first image
                axes = parsedheader["axes"]
                axesnames = [ax.name for ax in axes]
                self._add_grid_axis(axes[0], index=self.outimgdim[0])
                self._add_grid_axis(axes[1], index=self.outimgdim[1])

                # Stack axes: one value for each image
                stackaxes = {}
                stackaxisname = self.parameters["edfheader"]["stackvalue"]
                self._add_stack_axis(
                    stackaxisname,
                    parsedheader["stackvalue"].units,
                    index=self.outstackdim,
                )
                stackaxes[stackaxisname] = "stackvalue"
                for mot in self.parameters["edfheader"]["axesnamemap"]:
                    if mot not in axesnames and mot in parsedheader:
                        if not np.isnan(parsedheader[mot].magnitude):
                            stackaxes[mot] = mot
                            self._add_stack_axis(
                                mot, self.units.get(mot, "dimensionless")
                            )

                # Info axes: one or more values for each image
                self._add_info_axis("expotime", defaultunits=parsedheader["time"].units)
                infoaxes = {"expotime": "time"}
                if self.fluxnorm:
                    values = np.full(
                        (self.nstack, self.ndetfit), np.nan, dtype=np.float32
                    )
                    self._add_info_axis(
                        "xrfdetectorposition", values=values, defaultunits="cm"
                    )

                # Other scan info
                self.procinfo[axesnames[0]] = parsedheader.get("motslow", "")
                self.procinfo[axesnames[1]] = parsedheader.get("motfast", "")

            # Values for stack and info axes:
            # TODO: xrfdetectorposition
            for mot, hkey in stackaxes.items():
                self.axes[mot][imageindex] = parsedheader[hkey]
            for param, hkey in infoaxes.items():
                self.infoaxes[param][imageindex] = parsedheader[hkey]

            # Lazy add counters
            files = xiaimage.ctrfilenames_used(counters)
            files = xiaedf.xiagroupdetectors(files)

            n0, n1 = xiaimage.dshape[:2]

            def func(*args, **kw):
                img = nxlazy.readedf_func(*args, **kw)
                if img.shape != (n0, n1):
                    img2 = np.zeros((n0, n1), dtype=img.dtype)
                    m0, m1 = img.shape
                    m0 = min(m0, n0)
                    m1 = min(m1, n1)
                    img2[:m0, :m1] = img[:m0, :m1]
                    return img2
                else:
                    return img

            for detector, v1 in files.items():
                name = nxresult.Group(detector)

                # Prepare list of files
                if imageindex == 0:
                    self.stacks[name] = OrderedDict()
                    for ctr in v1:
                        self.stacks[name][ctr] = [None] * self.nstack

                for ctr, f in v1.items():
                    # Add counter file
                    o = nxlazy.LazyStackSlice(func=func)
                    o.appendarg_edf(f[0])
                    self.stacks[name][ctr][imageindex] = o
                    self.counters.add(ctr)

    def _add_grid_axis(self, axis, index=None):
        self.axes[axis.name] = axis
        if index is not None:
            self.axes_names[index] = axis.name

    def _add_stack_axis(self, name, unit, index=None):
        values = np.full(self.nstack, np.nan, dtype=np.float32)
        values = units.Quantity(values, units=unit)
        ax = axis.Axis(values, name=name)
        self._add_grid_axis(ax, index=index)

    def _add_info_axis(self, name, values=None, defaultunits=None, type=None):
        if values is None:
            values = np.full(self.nstack, np.nan, dtype=np.float32)
        u = self.units.get(name, defaultunits)
        values = units.Quantity(values, units=u)
        self.infoaxes[name] = axis.Axis(values, type=type, name=name)

    def _getscanparameters(self, header):
        """Get scan dimensions from header"""
        parameters = self.parameters
        o = spec.edfheader_parser(units=self.units, **parameters["edfheader"])
        metadata = parameters["metadata"]
        if metadata and metadata != "xia":
            defaultdims = None
        else:
            shape = self.xiastackraw.dshape
            defaultdims = shape[1], shape[2]
        return o.parse(header, defaultdims=defaultdims)

    def _stack_add_flux(self):
        if not self.fluxnorm or (
            "fluxcounter" not in self.parameters
            and "transmissioncounter" not in self.parameters
        ):
            return

        self._add_info_axis("i0_to_flux_offset", defaultunits="Hz")
        self._add_info_axis("i0_to_flux_factor", defaultunits="Hz")
        self._add_info_axis("it_to_flux_offset", defaultunits="Hz")
        self._add_info_axis("it_to_flux_factor", defaultunits="Hz")

        name = nxresult.Group(None)
        for imageindex in range(self.nstack):
            energy = self.axes[self.axes_names[self.outstackdim]][imageindex].magnitude
            time = self.infoaxes["refexpotime"][imageindex]

            if "fluxcounter" in self.parameters:
                op, _ = self.qxrfgeometry.I0op(
                    energy, expotime=time, removebeamfilters=False
                )
                if "calc_flux0" not in self.stacks[name]:
                    self.stacks[name]["calc_flux0"] = [None] * self.nstack
                o = nxlazy.LazyStackSlice(func=op)
                o.appendarg_h5dataset(
                    self.temp_nxresults[str(name)][self.parameters["fluxcounter"]]
                )
                self.stacks[name]["calc_flux0"][imageindex] = o
                self.infoaxes["i0_to_flux_offset"][imageindex] = op.b
                self.infoaxes["i0_to_flux_factor"][imageindex] = op.m

            if "transmissioncounter" in self.parameters:
                op, _ = self.qxrfgeometry.Itop(
                    energy, expotime=time, removebeamfilters=True
                )
                if "calc_fluxt" not in self.stacks[name]:
                    self.stacks[name]["calc_fluxt"] = [None] * self.nstack
                    self.stacks[name]["calc_transmission"] = [None] * self.nstack
                    self.stacks[name]["calc_absorbance"] = [None] * self.nstack

                o = nxlazy.LazyStackSlice(func=op)
                o.appendarg_h5dataset(
                    self.temp_nxresults[str(name)][
                        self.parameters["transmissioncounter"]
                    ]
                )
                self.stacks[name]["calc_fluxt"][imageindex] = o
                self.infoaxes["it_to_flux_offset"][imageindex] = op.b
                self.infoaxes["it_to_flux_factor"][imageindex] = op.m

                o = nxlazy.LazyStackSlice(func=nxlazy.transmission_func)
                o.appendarg_h5dataset(self.temp_nxresults[str(name)]["calc_fluxt"])
                o.appendarg_h5dataset(self.temp_nxresults[str(name)]["calc_flux0"])
                self.stacks[name]["calc_transmission"][imageindex] = o

                o = nxlazy.LazyStackSlice(func=nxlazy.absorbance_func)
                o.appendarg_h5dataset(
                    self.temp_nxresults[str(name)]["calc_transmission"]
                )
                self.stacks[name]["calc_absorbance"][imageindex] = o

    def _stack_add_xrffit(self):
        # Fit data and add elemental maps
        cfglist = self.parameters["pymcacfg"]
        ncfg = len(cfglist)
        if not ncfg:
            return

        logger.info("Fit XRF spectra ...")

        # Raises an error when the directory already exists
        outpath = self.outfitpath
        outpath.mkdir(force=False)
        if ncfg == 1:
            cfglist = cfglist * self.ndetfit
        else:
            cfglist = cfglist
            if len(cfglist) != self.ndetfit:
                raise RuntimeError(
                    "You need {} configuration files, {} provides.".format(
                        self.ndetfit, len(cfglist)
                    )
                )

        for imageindex, xiaimage in enumerate(self.xiastackproc):
            if self.fluxnorm:
                quantlist = [
                    {
                        "time": self.infoaxes["refexpotime"][imageindex]
                        .to("s")
                        .magnitude,
                        "flux": self.infoaxes["refflux"][imageindex].to("Hz").magnitude,
                        "area": self.infoaxes["activearea"][imageindex, i]
                        .to("cm**2")
                        .magnitude,
                        "anglein": self.infoaxes["anglein"][imageindex, i]
                        .to("deg")
                        .magnitude,
                        "angleout": self.infoaxes["angleout"][imageindex, i]
                        .to("deg")
                        .magnitude,
                        "distance": self.infoaxes["sampledetdistance"][imageindex, i]
                        .to("cm")
                        .magnitude,
                    }
                    for i in range(self.ndetfit)
                ]
            else:
                quantlist = [self.parameters["quantification"]] * self.ndetfit

            filestofit = xiaimage.datafilenames_used()
            filestofit = xiaedf.xiagroupdetectors(filestofit)
            filestofit = {nxresult.Group(k): v for k, v in filestofit.items()}
            for detector, cfg, quant in zip(self._detectors_to_fit, cfglist, quantlist):
                logger.info("Pymca fit of {} ...".format(detector))

                # Fit
                outname = "{}_{}_{:04d}_0000".format(
                    xiaimage.radix, detector.xialabel, xiaimage.mapnum
                )
                energy = (
                    self.axes[self.axes_names[self.outstackdim]][imageindex]
                    .to("keV")
                    .magnitude
                )
                if self.parameters["ignore_energy"]:
                    energy = None

                files, labels = PerformBatchFit(
                    filestofit[detector]["xia"],
                    outpath.path,
                    outname,
                    cfg,
                    energy,
                    fast=self.parameters["fastfitting"],
                    mlines=self.parameters["mlines"],
                    quant=quant,
                    addhigh=self.parameters["addhigh"],
                )

                # Prepare list of files
                name = nxresult.Group(detector)
                if imageindex == 0:
                    if name not in self.stacks:
                        self.stacks[name] = {}
                    for label in sorted(labels):
                        self.stacks[name][label] = [None] * self.nstack

                # Add file name
                for f, label in zip(files, labels):
                    o = nxlazy.LazyStackSlice()
                    o.appendarg_edf(f)
                    self.stacks[name][label][imageindex] = o

    def _postcorrect(self):
        if self._detectors_sumto:
            self._apply_postcorrect_xrf(self.fluxnormafter, self.dtcorafter)
            self._add_detectors()

            # if self.dtcorafter:
            #    self._apply_postcorrect_xrf(False, self.dtcorafter)
            # self._add_detectors()
            # if self.fluxnormafter:
            #    self._apply_postcorrect_xrf(self.fluxnormafter, False)
        else:
            self._apply_postcorrect_xrf(self.fluxnormafter, self.dtcorafter)

    def _add_detectors(self):
        # Datasets in sum groups that are already their shouldn't be changed:
        fixedgroups = {}

        # Add detector to sum and remove it:
        for k1, sumdest in self._detectors_sumto.items():
            if k1 not in self.stacks:
                continue

            if sumdest not in self.stacks:
                self.stacks[sumdest] = {}
            if sumdest not in fixedgroups:
                fixedgroups[sumdest] = list(self.stacks[sumdest].keys())

            logger.debug("Remove {} and add to {}".format(k1, sumdest))

            for k2 in self.stacks[k1]:
                if k2 in fixedgroups[sumdest]:
                    logger.debug("Do not add {}['{}'] to {}".format(k1, k2, sumdest))
                    continue
                if k2 not in self.stacks[sumdest]:
                    func = self._group_add_detectors(k2)
                    self.stacks[sumdest][k2] = [
                        nxlazy.LazyStackSlice(func=func, unpackargs=False)
                        for _ in range(self.nstack)
                    ]
                for imageindex, arg in enumerate(self.stacks[k1][k2]):
                    self.stacks[sumdest][k2][imageindex].appendarg(arg)

            self.stacks.pop(k1)

    def _group_add_detectors(self, grpname):
        if grpname in self.counters:
            func = nxlazy.sum_func
        elif grpname.startswith("w"):
            # mass fractions
            func = nxlazy.nanmean_func
        elif "chisq" in grpname:
            func = nxlazy.nanmax_func
        else:
            # peak areas
            func = nxlazy.sum_func

        return func

    def _apply_postcorrect_xrf(self, fluxnorm, dtcor):
        detectors = [detector for detector in self.stacks if detector.isdetector]
        normname = nxresult.Group(None)

        msg = "Corrections after XRF fitting:"
        if fluxnorm:
            msg += " flux normalization"
        if dtcor:
            if not msg.endswith(":"):
                msg += ","
            msg += " deadtime correction"
        if msg.endswith(":"):
            msg += " none"
        logger.info(msg)

        for k1 in detectors:
            for k2 in self.stacks[k1]:
                if k2 in self.counters:
                    continue

                keep = self.stacks[k1][k2]
                self.stacks[k1][k2] = [
                    nxlazy.LazyStackSlice(func=nxlazy.xrfnorm_func)
                    for _ in range(self.nstack)
                ]
                for imageindex, (arg, xiaimage) in enumerate(
                    zip(keep, self.xiastackproc)
                ):
                    # arguments: xrf,flux,fluxref,xiaimage
                    self.stacks[k1][k2][imageindex].appendarg(arg)
                    if fluxnorm:
                        self.stacks[k1][k2][imageindex].appendarg_h5dataset(
                            self.temp_nxresults[str(normname)]["calc_flux0"]
                        )
                        self.stacks[k1][k2][imageindex].appendarg(
                            self.infoaxes["refflux"][imageindex]
                        )
                    else:
                        self.stacks[k1][k2][imageindex].appendarg(None)
                        self.stacks[k1][k2][imageindex].appendarg(None)
                    if dtcor and not k1.issum:
                        self.stacks[k1][k2][imageindex].appendarg(xiaimage)
                        self.stacks[k1][k2][imageindex].appendarg(k1.number)
                    else:
                        self.stacks[k1][k2][imageindex].appendarg(None)
                        self.stacks[k1][k2][imageindex].appendarg(None)

    def _sort_stacks(self):
        # Sort stack on stack axis value
        mot = self.axes_names[self.outstackdim]
        ind = np.argsort(self.axes[mot].magnitude, kind="mergesort")
        self.axes[mot].values = self.axes[mot][ind]

        for name in self.infoaxes:
            self.infoaxes[name].values = self.infoaxes[name][ind, ...]

        for k1 in self.stacks:
            group = self.stacks[k1]
            for k2 in group:
                group[k2] = [group[k2][i] for i in ind]

    def _exportgroups(self):
        """Export groups of EDF stacks, summed or not"""
        outshape = None
        for k1 in self.stacks:  # detector or counter group
            nxdata = self.temp_nxresults.nxdata(str(k1))
            # stack subgroup (Al-K, S-K, xmap_x1c, ...)
            for k2 in self.stacks[k1]:
                signal = nxdata[str(k2)]
                for imageindex, slicedata in enumerate(self.stacks[k1][k2]):
                    if not signal.exists:
                        logger.info("saving {}/{}".format(k1, k2))
                    logger.debug(
                        "Slice generator (index {}): {}".format(imageindex, slicedata)
                    )

                    data = slicedata.data(imageindex, self.outstackdim)

                    # Allocate destination
                    if not signal.exists:
                        if outshape is None:
                            outshape = [0, 0, 0]
                            outshape[self.outimgdim[0]] = data.shape[0]
                            outshape[self.outimgdim[1]] = data.shape[1]
                            outshape[self.outstackdim] = len(self.stacks[k1][k2])
                        nxdata.add_signal(
                            name=signal.name,
                            shape=outshape,
                            chunks=True,
                            dtype=np.float32,
                        )

                    # Some rows too much or rows missing:
                    if outshape[self.outimgdim[0]] > data.shape[0]:
                        data = np.pad(
                            data,
                            ((0, outshape[self.outimgdim[0]] - data.shape[0]), (0, 0)),
                            "constant",
                            constant_values=0,
                        )
                    elif outshape[self.outimgdim[0]] < data.shape[0]:
                        data = data[0 : outshape[self.outimgdim[0]], :]

                    # Save data
                    with signal.open() as dset:
                        if self.outstackdim == 0:
                            dset[imageindex, ...] = data
                        elif self.outstackdim == 1:
                            dset[:, imageindex, :] = data
                        else:
                            dset[..., imageindex] = data

                if signal.exists:
                    self.stacks[k1][k2] = signal

            if list(nxdata.signals):
                nxdata.set_axes(*self.axes_names)

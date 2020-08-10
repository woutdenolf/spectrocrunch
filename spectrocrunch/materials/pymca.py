# -*- coding: utf-8 -*-

from __future__ import absolute_import

import contextlib
import re
import copy
import operator
from types import MethodType
import collections

from . import mixture
from . import element
from . import compoundfromformula
from . import xrayspectrum
from . import compoundfromlist
from . import types
from ..utils import units
from ..utils import instance
from ..patch import xraylib
from ..io.mca import save as savemca
from ..utils.copyable import Copyable

import numpy as np
import scipy.interpolate
import scipy.optimize
import matplotlib.pyplot as plt
import fisx
from PyMca5.PyMcaPhysics.xrf import ClassMcaTheory
from PyMca5.PyMcaPhysics.xrf import ConcentrationsTool
from PyMca5.PyMcaIO import ConfigDict

try:
    from silx.gui import qt
    from PyMca5.PyMcaGui.physics.xrf import McaAdvancedFit
except ImportError:
    qt = None
    McaAdvancedFit = None


class PymcaBaseHandle(Copyable):
    def __init__(self):
        self.mcafit = ClassMcaTheory.McaTheory()
        self.ctool = ConcentrationsTool.ConcentrationsTool()
        self.app = None

    def loadfrompymca(self, filename=None, config=None):
        if filename is not None:
            fconfig = ConfigDict.ConfigDict()
            fconfig.read(filename)
            self.mcafit.configure(fconfig)
        if config is not None:
            self.mcafit.configure(config)

    def savepymca(self, filename):
        self.mcafit.config.write(filename)

    def configfromfitresult(self, digestedresult):
        config = copy.deepcopy(digestedresult["config"])
        # Global parameters
        parammap = {
            "Zero": ("detector", "zero"),
            "Gain": ("detector", "gain"),
            "Noise": ("detector", "noise"),
            "Fano": ("detector", "fano"),
            "Sum": ("detector", "sum"),
            "ST AreaR": ("peakshape", "st_arearatio"),
            "ST SlopeR": ("peakshape", "st_sloperatio"),
            "LT AreaR": ("peakshape", "lt_arearatio"),
            "LT SlopeR": ("peakshape", "lt_sloperatio"),
            "STEP HeightR": ("peakshape", "step_heightratio"),
            "Eta Factor": ("peakshape", "eta_factor"),
        }
        nparams = len(digestedresult["parameters"])
        ngroups = len(digestedresult["groups"])
        nglobal = nparams - ngroups
        for i in range(nglobal):
            param = digestedresult["parameters"][i]
            if param in parammap:
                key = parammap[param]
                config[key[0]][key[1]] = digestedresult["fittedpar"][i]
        return config

    @staticmethod
    def _parse_fitresult_massfractions(wfrac1, wfrac0, exclude):
        # Group per element
        wfrac2 = {}
        for zgroup, w in wfrac1.items():
            element = zgroup.element
            if element in exclude:
                continue
            group = zgroup.element.shells[0]
            if element not in wfrac2:
                wfrac2[element] = {}
            wfrac2[element][group] = w
        # Select mass fraction from most reliable group
        wfrac3 = {}
        for element, groupdict in wfrac2.items():
            w = sorted(groupdict.items(), key=operator.itemgetter(0))[0][1]
            wfrac3[element] = w
        # Modify original mass fractions
        wfracr = dict(wfrac0)
        snew = 1 - sum(w for element, w in wfrac0.items() if element not in wfrac3)
        sold = sum(wfrac3.values())
        norm = snew / sold
        for k, v in wfrac3.items():
            if v:
                wfracr[k] = v * norm
        return wfracr

    def _pymcaconfig_massfractions(self, config):
        """
        Extract elemental mass fractions from the pymca config sample
        """
        materials = collections.OrderedDict()
        enabled, name, density = config["attenuators"]["Matrix"][:3]
        if name == "MULTILAYER":
            for i in range(len(config["multilayer"])):
                enabled, name, density = config["multilayer"]["Layer{}".format(i)][:3]
                if enabled:
                    mat = self.loadfrompymca_material(config, name, density)
                    materials[name] = mat.elemental_massfractions()
        elif enabled:
            mat = self.loadfrompymca_material(config, name, density)
            materials[name] = mat.elemental_massfractions()
        return materials

    def configadjustconcentrations(self, config, lmassfractions, exclude=None):
        """
        Modify pymca config sample materials based on list of elemental concentrations
        """
        if not exclude:
            exclude = []
        materials = self._pymcaconfig_massfractions(config)
        if not materials:
            return
        if len(materials) != len(lmassfractions):
            raise RuntimeError("Number of layers in config and fitresult do not match")
        for (name, wfrac0), wfrac1 in zip(materials.items(), lmassfractions):
            wfrac2 = self._parse_fitresult_massfractions(wfrac1, wfrac0, exclude)
            CompoundList = [str(k) + "1" for k in wfrac2.keys()]
            CompoundFraction = list(wfrac2.values())
            config["materials"][name]["CompoundList"] = CompoundList
            config["materials"][name]["CompoundFraction"] = CompoundFraction

    def setdata(self, y):
        x = np.arange(len(y), dtype=np.float32)
        self.mcafit.setData(x, y)

    def fit(
        self,
        loadfromfit=False,
        concentrations=True,
        spectra=True,
        environ_elements=None,
    ):
        # Fit
        self.mcafit.estimate()
        fitresult, digestedresult = self.mcafit.startfit(digest=1)
        # Background
        yback = digestedresult.pop("continuum")  # polynomial + snip
        if self.mcafit.STRIP:
            ysnip = np.ravel(self.mcafit.zz)
        else:
            ysnip = 0
        digestedresult["yback"] = yback
        digestedresult["ysnip"] = ysnip
        # Parse result
        result = self.miscfromfitresult(digestedresult)
        if concentrations:
            self.concentrationsfromfitresult(digestedresult, out=result)
        if spectra:
            self.spectrafromfitresult(digestedresult, out=result)
        # Load parameters from fit (if you did "load from fit")
        if loadfromfit:
            config = self.configfromfitresult(digestedresult)
            if concentrations:
                self.configadjustconcentrations(
                    config, result["lmassfractions"], exclude=environ_elements
                )
            self.loadfrompymca(config=config)
        return result

    def fitgui(
        self,
        ylog=False,
        legend="data",
        loadfromfit=False,
        concentrations=True,
        spectra=True,
        environ_elements=None,
    ):
        if self.app is None:
            self.app = qt.QApplication([])
        w = McaAdvancedFit.McaAdvancedFit()

        # Copy mcafit
        x = self.mcafit.xdata0
        y = self.mcafit.ydata0
        w.setData(x, y, legend=legend, xmin=0, xmax=len(y) - 1)
        w.mcafit.configure(self.mcafit.getConfiguration())

        # GUI for fitting
        if ylog:
            w.graphWindow.yLogButton.click()
        w.graphWindow.energyButton.click()
        # w.graphWindow.setGraphYLimits(min(y[y>0]),None)
        w.refreshWidgets()
        w.show()
        result = self.app.exec_()

        # Parse result
        result = {}
        try:
            digestedresult = w.mcafit.digestresult()
        except AttributeError:
            digestedresult = None

        if digestedresult:
            # Load parameters from fit (if you did "load from fit")
            self.miscfromfitresult(digestedresult, out=result)
            if concentrations:
                self.concentrationsfromfitresult(digestedresult, out=result)
            if spectra:
                self.spectrafromfitresult(digestedresult, out=result)
            if loadfromfit:
                config = w.mcafit.getConfiguration()
                if concentrations:
                    self.configadjustconcentrations(
                        config, result["lmassfractions"], exclude=environ_elements
                    )
                self.loadfrompymca(config=config)
        else:
            if concentrations or spectra:
                result = self.fit(
                    loadfromfit=loadfromfit,
                    concentrations=concentrations,
                    spectra=spectra,
                )

        return result

    def processfitresult(self, digestedresult, originalconcentrations=False):
        ctoolcfg = self.ctool.configure()
        ctoolcfg.update(digestedresult["config"]["concentrations"])
        return self.ctool.processFitResult(
            config=ctoolcfg,
            fitresult={"result": digestedresult},
            elementsfrommatrix=originalconcentrations,
            fluorates=self.mcafit._fluoRates,
            addinfo=True,
        )

    def concentrationsfromfitresult(self, digestedresult, out=None):
        # Mass fractions:
        #
        # For a group (e.g Fe-K) and i loops over the layers:
        #  grouparea = flux.time.solidangle/(4.pi).sum_i[massfrac_i.grouprate_i]
        #
        # When Fe is present in only one layer j:
        #  massfrac_j = grouparea/(flux.time.solidangle/(4.pi).grouprate_i)
        #
        # When Fe present in multiple layers:
        #  grouparea = flux.time.solidangle/(4.pi).massfrac_avg.sum_i[b_i.grouprate_i]
        # where b_i is 1 or 0 depending on whether it is present in the matrix definition
        # When the element is not in the matrix definition b_i=1 for all layers.
        #
        #  massfrac_avg = grouparea/(flux.time.solidangle/(4.pi).sum_i[b_i.grouprate_i])
        if out is None:
            out = {}
        conresult, addinfo = self.processfitresult(
            digestedresult, originalconcentrations=False
        )

        out["massfractions"] = {
            element.Element.fluozgroup(k): v
            for k, v in conresult["mass fraction"].items()
        }

        if "layerlist" in conresult:
            nlayers = len(conresult["layerlist"])
            if nlayers > 0:
                out["lmassfractions"] = [
                    {
                        element.Element.fluozgroup(k): v
                        for k, v in conresult[k]["mass fraction"].items()
                    }
                    for k in conresult["layerlist"]
                ]
            else:
                nlayers = 1
                out["lmassfractions"] = []
        else:
            nlayers = 1
            out["lmassfractions"] = []

        if len(out["lmassfractions"]) == 0:
            out["lmassfractions"] = [out["massfractions"]]

        # Group rates:
        #   grouprate = solidangle/(4.pi).sum_i[massfrac_i.grouprate_i]   (i loops over the layers)
        out["rates"] = {}  # Ifluo / (Flux.time)
        safrac = addinfo["SolidAngle"]
        # assert(self.sample.geometry.solidangle/(4*np.pi)==addinfo["SolidAngle"])
        # assert(self.flux*self.time==addinfo["I0"])
        for zgroup in out["massfractions"]:
            ele = str(zgroup.element)
            group = str(zgroup.group)
            grouprate = 0.0
            for layer in range(1, nlayers + 1):
                if ele in self.mcafit._fluoRates[layer]:
                    massfrac_l = self.mcafit._fluoRates[layer][ele]["mass fraction"]
                    grouprate_l = self.mcafit._fluoRates[layer][ele]["rates"][
                        "{} xrays".format(group)
                    ]
                    grouprate += massfrac_l * grouprate_l
            grouprate *= safrac
            out["rates"][zgroup] = grouprate

        # Fitted areas (elements)
        out["fitareas"] = {
            element.Element.fluozgroup(k): v for k, v in conresult["fitarea"].items()
        }

        # Fitted areas (scatter)
        scatterpeaks = {
            "Compton": {"energy": [], "fitarea": []},
            "Rayleigh": {"energy": [], "fitarea": []},
        }
        nen = len(instance.asarray(digestedresult["energy"]))
        for k, v in digestedresult.items():
            if k.startswith("Scatter"):
                if "Compton" in k:
                    peak = "Compton"
                else:
                    peak = "Rayleigh"
                energy, fitarea = zip(
                    *[(v[k2]["energy"], v[k2]["fitarea"]) for k2 in v["peaks"]]
                )
                scatterpeaks[peak]["energy"] += energy
                scatterpeaks[peak]["fitarea"] += fitarea
        for peak in scatterpeaks:
            energy = np.array(scatterpeaks[peak]["energy"])
            fitarea = np.array(scatterpeaks[peak]["fitarea"])
            if peak == "Compton":
                k = xrayspectrum.ComptonLine(energy)
            else:
                k = xrayspectrum.RayleighLine(energy)
            out["fitareas"][k] = fitarea

        # Fitted rates
        out["fitrates"] = {k: v / addinfo["I0"] for k, v in out["fitareas"].items()}

        # self._print_pymcainternals_rates()
        return out

    def _pymcainternals_solidanglefrac(self):
        radius2 = self.mcafit.config["concentrations"]["area"] / np.pi
        distance = self.mcafit.config["concentrations"]["distance"]
        return 0.5 * (1.0 - (distance / np.sqrt(distance ** 2 + radius2)))

    def _print_pymcainternals_rates(self):
        safrac = self._pymcainternals_solidanglefrac()

        rowfmt = "{:>6}{:>8}{:>20}{:>10}{:>10}{:>20}"
        print(rowfmt.format("Layer", "Element", "MassFrac", "Line", "Energy", "Rate"))
        for i, layer in enumerate(self.mcafit._fluoRates):
            if i == 0:
                continue
            for ele in sorted(layer):
                lines = []
                for k, v in layer[ele].items():
                    if k.endswith("xrays"):
                        lines.extend(v)
                for line in lines:
                    w = layer[ele]["mass fraction"]
                    rate = layer[ele][line]["rate"] * w * safrac
                    energy = layer[ele][line]["energy"]
                    if rate > 0:
                        print(rowfmt.format(i - 1, ele, w, line, energy, rate))

    def miscfromfitresult(self, digestedresult, out=None):
        if out is None:
            out = {}
        out["chisq"] = digestedresult["chisq"]
        return out

    def spectrafromfitresult(self, digestedresult, out=None):
        conresult, addinfo = self.processfitresult(
            digestedresult, originalconcentrations=True
        )

        # Group rates:
        if False:
            # TODO: set rates of elements not in matrix to zero
            out["rates"] = {}
            for group in out["area"]:
                out["rates"][element.fluozgroup(group)] = (
                    conresult["area"][group] / addinfo["I0"]
                )

        # Matrix spectrum
        nparams = len(digestedresult["parameters"])
        ngroups = len(digestedresult["groups"])
        nglobal = nparams - ngroups
        parameters = list(digestedresult["fittedpar"])
        for i, group in enumerate(digestedresult["groups"]):
            if group in conresult["area"]:
                # Replace fitted with theoretical peak area
                parameters[nglobal + i] = conresult["area"][group]
            else:
                # Scattering peaks don't have a theoretical peak area
                parameters[nglobal + i] = 0.0
        ymatrix = self.mcafit.mcatheory(parameters, digestedresult["xdata"])

        yback = digestedresult["yback"]  # polynomial + snip
        ysnip = digestedresult["ysnip"]

        def interpol(x, spectrum):
            return scipy.interpolate.interp1d(
                x,
                spectrum,
                kind="nearest",
                bounds_error=False,
                fill_value=(spectrum[0], spectrum[-1]),
            )

        def interpol_energy(prof):
            return interpol(digestedresult["energy"], prof)

        def interpol_channel(prof):
            return interpol(digestedresult["xdata"], prof)

        if out is None:
            out = {}
        out["energy"] = digestedresult["energy"]
        out["channels"] = digestedresult["xdata"]
        out["y"] = digestedresult["ydata"]
        out["yfit"] = digestedresult["yfit"]
        out["yback"] = yback
        out["interpol_energy"] = interpol_energy
        out["interpol_channel"] = interpol_channel
        out["ypileup"] = digestedresult["pileup"]
        # polynomial is already in ymatrix, snip not
        out["ymatrix"] = ymatrix + ysnip

        def _plot(ylog=False, label="data"):
            plt.plot(out["energy"], out["y"], "+", label=label)
            plt.plot(out["energy"], out["yfit"], label=label + " fit")
            plt.plot(out["energy"], out["ymatrix"], label=label + " matrix")
            plt.plot(out["energy"], out["yback"], label=label + " background")
            if ylog:
                plt.gca().set_yscale("log", basey=10)
            plt.legend()
            plt.xlabel("MCA (keV)")
            plt.ylabel("Counts")

        out["plot"] = _plot

        return out


class PymcaHandle(PymcaBaseHandle):
    """Class that converts spectrocrunch to pymca parameters
       and vice versa. Also used to call the pymca fitting
       API for single spectra.
    """

    def __init__(
        self,
        sample=None,
        emin=None,
        emax=None,
        energy=None,
        weights=None,
        scatter=None,
        flux=1e9,
        time=0.1,
        escape=True,
        pileup=False,
        ninteractions=1,
        linear=False,
        snip=True,
        continuum=0,
        SingleLayerStrategy=None,
        noisepropagation=True,
    ):
        self.sample = sample
        self.set_source(energy=energy, weights=weights, scatter=scatter)
        self.emin = emin
        self.emax = emax
        self.linear = linear
        self.escape = escape
        self.pileup = pileup
        self.snip = snip
        self.continuum = continuum
        self.noisepropagation = noisepropagation
        self.ninteractions = ninteractions
        self.SingleLayerStrategy = SingleLayerStrategy
        self.flux = flux
        self.time = time
        super(PymcaHandle, self).__init__()

    def __str__(self):
        s = zip(
            instance.asarray(self.energy),
            instance.asarray(self.weights),
            instance.asarray(self.scatter),
        )
        s = "\n ".join(
            "{} keV: {} % (Scatter: {})".format(k, v * 100, sc) for k, v, sc in s
        )
        if self.sample:
            sample = str(self.sample)
            geometry = str(self.sample.geometry)
        else:
            sample = geometry = None
        return "Flux = {:~e}\nTime = {:~}\nSource lines:\n {}\n{}\n{}".format(
            self.flux, self.time, s, sample, geometry
        )

    def set_source(self, energy=None, weights=None, scatter=None):
        if energy is None:
            return
        self.energy = units.umagnitude(energy, "keV")
        if weights is None:
            self.weights = np.ones_like(self.energy)
        else:
            self.weights = weights
        if scatter is None:
            self.scatter = np.ones_like(self.energy)
        elif isinstance(scatter, bool):
            if scatter:
                self.scatter = np.ones_like(self.energy)
            else:
                self.scatter = np.zeros_like(self.energy)
        else:
            self.scatter = scatter

    @property
    def flux(self):
        return self._flux

    @flux.setter
    def flux(self, value):
        self._flux = units.Quantity(value, "Hz").to("Hz")

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, value):
        self._time = units.Quantity(value, "s").to("s")

    @property
    def I0(self):
        return (self.flux * self.time).to("dimensionless").magnitude

    @property
    def emax_strict(self):
        return np.max(self.energy)

    @property
    def emax(self):
        if self._emax is None:
            # include elastic scattering peaks
            p = np.max(self.energy)
            if self.sample:
                s = np.sqrt(self.sample.geometry.detector.gaussianVAR(p))
            else:
                s = 0.1
            return p + 3 * s
        else:
            return self._emax

    @emax.setter
    def emax(self, value):
        self._emax = value

    @property
    def emin(self):
        if self._emin is None:
            # include Al-K lines
            p = xraylib.LineEnergy(13, xraylib.KL3_LINE)
            if self.sample:
                s = np.sqrt(self.sample.geometry.detector.gaussianVAR(p))
            else:
                s = 0.1
            return p - 3 * s
        else:
            return self._emin

    @emin.setter
    def emin(self, value):
        self._emin = value

    def addtopymca_lstsq(self, cfg):
        cfg["fit"]["maxiter"] = 500

    def addtopymca_beam(self, cfg):
        # TODO: move to source?
        cfg["fit"]["energy"] = instance.asarray(self.energy).tolist()
        cfg["fit"]["energyweight"] = instance.asarray(self.weights).tolist()
        cfg["fit"]["energyscatter"] = (
            instance.asarray(self.scatter).astype(int).tolist()
        )
        cfg["fit"]["energyflag"] = np.ones_like(
            cfg["fit"]["energy"], dtype=int
        ).tolist()
        cfg["fit"]["scatterflag"] = int(any(cfg["fit"]["energyscatter"]))

        # Not just the source:
        cfg["concentrations"]["flux"] = self.flux.to("Hz").magnitude
        cfg["concentrations"]["time"] = self.time.to("s").magnitude

    def loadfrompymca_beam(self, cfg):
        ind = np.asarray(cfg["fit"]["energyflag"], dtype=bool).tolist()

        self.energy = instance.asarray(cfg["fit"]["energy"])
        self.energy = self.energy[ind].astype(float)
        self.weights = np.asarray(cfg["fit"]["energyweight"], dtype=float)
        self.weights = self.weights[ind]
        self.scatter = np.asarray(cfg["fit"]["energyscatter"], dtype=int)
        self.scatter = self.scatter[ind]
        self.flux = cfg["concentrations"]["flux"]
        self.time = cfg["concentrations"]["time"]

    def addtopymca_background(self, cfg):
        # bmodelbkg = self.sample.geometry.detector.bstail or \
        #            self.sample.geometry.detector.bltail or \
        #            self.sample.geometry.detector.bstep
        cfg["fit"]["stripflag"] = int(self.snip)
        cfg["fit"]["stripalgorithm"] = 1
        cfg["fit"]["snipwidth"] = 100
        cfg["fit"]["continuum"] = int(self.continuum)

    def loadfrompymca_background(self, cfg):
        self.snip = bool(cfg["fit"]["stripflag"])
        self.continuum = int(cfg["fit"]["continuum"])

    def addtopymca_other(self, cfg):
        cfg["fit"]["escapeflag"] = int(self.escape)
        cfg["fit"]["sumflag"] = int(self.pileup)
        cfg["fit"]["linearfitflag"] = int(self.linear)
        cfg["concentrations"]["usemultilayersecondary"] = self.ninteractions - 1
        cfg["fit"]["fitweight"] = int(self.noisepropagation)

    def addtopymca_custom(self, cfg):
        pass

    def addcustomnpymcamethod(self, method):
        self.addtopymca_custom = MethodType(method, self)

    def addtopymca_strategy(self, cfg):
        strategy = "SingleLayerStrategy"
        dic = getattr(self, strategy)
        if dic:
            dic = dict(dic)
            cfg["fit"]["strategy"] = strategy
            cfg["fit"]["strategyflag"] = 1
            if strategy not in cfg:
                cfg[strategy] = {}
                cfg[strategy]["layer"] = "Auto"
                cfg[strategy]["iterations"] = 3

            def func(mat):
                return mat if instance.isstring(mat) else mat.pymcaname

            cfg[strategy]["completer"] = func(dic.pop("completer", "-"))
            dic2 = {}
            for el, mat in dic.items():
                if el in cfg["peaks"]:
                    k = "{} {}".format(el, sorted(cfg["peaks"][el])[0])
                    dic2[k] = func(mat)
            cfg[strategy]["peaks"] = dic2.keys()
            cfg[strategy]["materials"] = dic2.values()
            cfg[strategy]["flags"] = [1] * len(dic2)
        else:
            cfg["fit"]["strategyflag"] = 0

    def loadfrompymca_other(self, cfg):
        self.escape = bool(cfg["fit"]["escapeflag"])
        self.pileup = bool(cfg["fit"]["sumflag"])
        self.linear = bool(cfg["fit"]["linearfitflag"])
        self.ninteractions = cfg["concentrations"]["usemultilayersecondary"] + 1
        self.noisepropagation = bool(cfg["fit"]["fitweight"])

    def addtopymca_material(self, cfg, material, defaultthickness=1e-4):
        return material.topymca(cfg, defaultthickness=defaultthickness)

    def loadfrompymca_material(self, cfg, matname, density):
        formula = re.compile("^(([A-Z][a-z]?)([0-9]*))+$")
        purelement = re.compile("^(?P<element>[A-Z][a-z]?)1?$")

        # Material name is a formula:
        name = mixture.Mixture.namefrompymca(matname)
        mf = formula.match(name)
        if mf:
            me = purelement.match(name)
            if me:
                return element.Element(me.groupdict()["element"])
            else:
                return compoundfromformula.CompoundFromFormula(name, density)

        # Parse material info:
        dic = cfg["materials"][matname]
        if not instance.isarray(dic["CompoundList"]):
            dic["CompoundList"] = [dic["CompoundList"]]
        if not instance.isarray(dic["CompoundFraction"]):
            dic["CompoundFraction"] = [dic["CompoundFraction"]]
        massfractions = np.array(dic["CompoundFraction"]) / sum(dic["CompoundFraction"])
        density = dic["Density"]

        # List of materials:
        lst = [
            self.loadfrompymca_material(cfg, mat, density)
            for mat in dic["CompoundList"]
        ]

        # All elements:
        if all([isinstance(e, element.Element) for e in lst]):
            return compoundfromlist.CompoundFromList(
                lst, massfractions, types.fraction.mass, density, name=name
            )
        else:
            return mixture.Mixture(lst, massfractions, types.fraction.mass, name=name)

    def loadfrompymca(self, filename=None, config=None):
        """Update pymca parameters and convert to spectrocrunch
        """
        super(PymcaHandle, self).loadfrompymca(filename=filename, config=config)
        config = self.mcafit.getConfiguration()
        self.loadfrompymca_beam(config)
        self.loadfrompymca_background(config)
        self.loadfrompymca_other(config)
        self.sample.loadfrompymca(self, config)

    def addtopymca(self, fresh=False, onlygeometry=False):
        """Convert spectrocrunch to pymca parameters
        """
        if fresh:
            config = self.mcafit.getStartingConfiguration()
            config["attenuators"] = {}
            try:
                self.mcafit.useFisxEscape(True)
            except:
                pass
            self.mcafit.disableOptimizedLinearFit()
        else:
            config = self.mcafit.getConfiguration()

        if onlygeometry:
            self.sample.geometry.addtopymca(self, config)
        else:
            self.addtopymca_beam(config)
            self.addtopymca_background(config)
            self.addtopymca_lstsq(config)
            self.addtopymca_other(config)
            self.sample.addtopymca(self, config)
            self.addtopymca_strategy(config)
        self.addtopymca_custom(config)
        self.mcafit.configure(config)

    def configurepymca(self, **kwargs):
        """Convert spectrocrunch to pymca parameters
        """
        self.addtopymca(**kwargs)

    def savepymca(self, filename):
        """Convert spectrocrunch to pymca parameters and then save
        """
        self.configurepymca()
        super(PymcaHandle, self).savepymca(filename)

    @contextlib.contextmanager
    def fitenergy_context(self):
        # Keep original config
        self.configurepymca()
        config = self.mcafit.getConfiguration()
        keep = self.emin, self.emax, self.energy

        # Change config (extend emin and emax)
        emin = np.min(self.energy)
        line = xrayspectrum.ComptonLine(emin)
        emin = line.energy(polar=np.pi)
        s = np.sqrt(self.sample.geometry.detector.gaussianVAR(emin))
        emin = min(keep[0], emin - 5 * s)
        emin = max(emin, 0)

        emax = np.max(self.energy)
        s = np.sqrt(self.sample.geometry.detector.gaussianVAR(emax))
        emax = max(keep[1], emax + 3 * s)

        self.emin, self.emax = emin, emax
        self.configurepymca()

        yield

        # Reset config
        self.emin, self.emax, self.energy = keep
        self.mcafit.configure(config)
        self.configurepymca()

    def fitenergy(self, plot=False):
        with self.fitenergy_context():
            # Fit reduced range
            result = self.fit()
            if plot:
                plt.plot(result["energy"], result["y"], label="data")
                plt.plot(result["energy"], result["yfit"], label="pymca fit")
                plt.legend()
                plt.show()

            config = self.mcafit.getConfiguration()

            # Repeat fit
            def func(params):
                params = instance.asarray(params)
                self.energy = params[:-1]
                self.sample.geometry.angleout = params[-1]
                self.configurepymca()

                self.mcafit.configure(config)

                fitresult, digestedresult = self.mcafit.startfit(digest=1)
                return digestedresult["yfit"] - digestedresult["ydata"]

            guess = instance.asarray(self.energy).tolist()
            guess.append(self.sample.geometry.angleout)

            params, success = scipy.optimize.leastsq(func, guess)

        success = success > 0 and success < 5
        if success:
            params = instance.asarray(params)
            self.energy = params[:-1]
            self.sample.geometry.angleout = params[-1]
            self.configurepymca()

    def xrayspectrum(self, **kwargs):
        return self.sample.xrayspectrum(
            self.energy,
            emin=self.emin,
            emax=self.emax,
            weights=self.weights,
            ninteractions=self.ninteractions,
            **kwargs
        )

    def xraygrouprates(self, **kwargs):
        spectrum = self.xrayspectrum(**kwargs)
        groups = spectrum.linegroups(convert=True)
        pattern = "^(?P<element>[A-Z][a-z]?)-(?P<group>[A-Z]).?$"
        groups2 = {}
        for k in groups:
            m = re.match(pattern, k)
            if m:
                m = m.groupdict()
                g = "{}-{}".format(m["element"], m["group"])
                A = sum(groups[k].values())
                if g in groups2:
                    groups2[g] += A
                else:
                    groups2[g] = A
            else:
                groups2[k] = sum(groups[k].values())

        return groups2

    def xraygroupareas(self, **kwargs):
        return {k: v * self.I0 for k, v in self.xraygrouprates(**kwargs).items()}

    def mca(self, histogram=True, decomposed=0, energy=False, full=False, **kwargs):
        spectrum = self.xrayspectrum(**kwargs)
        if decomposed == 0:
            x, y, ylabel = spectrum.sumspectrum(fluxtime=self.I0, histogram=histogram)
        elif decomposed == 1:
            x, y, ylabel = spectrum.peakspectrum(
                fluxtime=self.I0, histogram=histogram, group=True
            )
        else:
            x, y, ylabel = spectrum.peakspectrum(
                fluxtime=self.I0, histogram=histogram, group=False
            )
        if not energy:
            a, b = spectrum.channellimits
            nchan = int(2 ** np.ceil(np.log2(b + 1)))
            mca = np.zeros(nchan, dtype=y.dtype)
            mca[a : b + 1] = y
            y = mca
            x = np.arange(nchan)
        if full:
            return x, y, ylabel
        else:
            return y

    def savemca(self, filename, mode="w", func=None, **kwargs):
        y = self.mca(**kwargs)
        if instance.iscallable(func):
            y = func(y)
        savemca(
            y,
            filename,
            mode=mode,
            zero=self.sample.geometry.detector.mcazero,
            gain=self.sample.geometry.detector.mcagain,
        )


class FisxConfig:

    FISXMATERIALS = None

    def init(self):
        if self.FISXMATERIALS is None:
            self.FISXMATERIALS = fisx.Elements()
            self.FISXMATERIALS.initializeAsPyMca()

    @contextlib.contextmanager
    def init_ctx(self):
        self.init()
        yield

    def setDetector(self, setup, detector):
        with self.init_ctx():
            detector.addtofisx(setup, self.FISXMATERIALS)

    def addtofisx_material(self, material, defaultthickness=1e-4):
        with self.init_ctx():
            return material.tofisx(
                self.FISXMATERIALS, defaultthickness=defaultthickness
            )

# -*- coding: utf-8 -*-

import PyMca5
from PyMca5.PyMcaPhysics.xrf import McaAdvancedFitBatch
from PyMca5.PyMcaPhysics.xrf import FastXRFLinearFit
from PyMca5.PyMcaPhysics.xrf import ClassMcaTheory
from PyMca5.PyMca import EDFStack
from PyMca5.PyMcaIO import ConfigDict

try:
    from PyMca5.PyMcaPhysics.xrf.McaAdvancedFitBatch import (
        OutputBuffer as OutputBufferBase,
    )
except ImportError:
    OutputBuffer = None
else:

    class OutputBuffer(OutputBufferBase):
        @property
        def outputDirLegacy(self):
            return self.outputDir


import numpy as np
import re
import os
import glob
from contextlib import contextmanager
import matplotlib.pyplot as plt
import logging

from ..utils import instance
from ..io import edf
from ..io import localfs
from ..io import utils as ioutils

logger = logging.getLogger(__name__)


def AdaptPyMcaConfigFile(filename, *args, **kwargs):
    cfg = ConfigDict.ConfigDict(filelist=filename)
    AdaptPyMcaConfig(cfg, *args, **kwargs)
    cfg.write(filename)


@contextmanager
def tempPyMcaConfigFile(cfg):
    ioutils.temporary_filename
    with ioutils.TemporaryFilename(suffix=".cfg") as filename:
        cfg.write(filename)
        yield filename


def AdaptPyMcaConfig_energy(cfg, energy, addhigh):
    if energy is None or not np.isfinite(energy):
        return

    # Extract source lines
    ind = instance.asarray(cfg["fit"]["energyflag"]).astype(bool)
    norg = len(ind)
    nenergies = ind.sum() + bool(addhigh)

    def extract(name, default=np.nan):
        arr = cfg["fit"][name]
        if instance.isarray(arr):
            arr = [instance.asnumber(v) for v in arr]
        arr = instance.asarray(arr)

        # Select based on energyflag
        narr = len(arr)
        if narr < norg:
            arr = np.append(arr, [default] * (norg - narr))
        arr = arr[0:norg][ind]

        # At least nenergies
        narr = len(arr)
        if narr < nenergies:
            arr = np.append(arr, [default] * (nenergies - narr))
        return arr

    cfg_energy = extract("energy", default=np.nan)
    cfg_energyweight = extract("energyweight", default=np.nan)
    cfg_energyflag = extract("energyflag", default=1)
    cfg_energyscatter = extract("energyscatter", default=0)

    # Modify energy
    cfg_energy = cfg_energy / cfg_energy[0] * energy
    cfg_energyweight = cfg_energyweight / cfg_energyweight[0]

    # Add missing lines
    for i in range(nenergies):
        if not np.isfinite(cfg_energy[i]):
            if i == 0:
                cfg_energy[i] = energy
            else:
                cfg_energy[i] = addhigh * energy
        if not np.isfinite(cfg_energyweight[i]):
            if i == 0:
                cfg_energyweight[i] = 1
            else:
                cfg_energyweight[i] = 1e-10

    # Remove extract line when it was already there
    if addhigh:
        if (
            cfg_energyweight[-2] / cfg_energyweight[0] < 1e-5
            and cfg_energy[-2] > energy
        ):
            nenergies -= 1
            cfg_energy = cfg_energy[:-1]
            cfg_energyweight = cfg_energyweight[:-1]
            cfg_energyflag = cfg_energyflag[:-1]
            cfg_energyscatter = cfg_energyscatter[:-1]

    # List with original size
    def reset(arr, default=0):
        arr = arr.tolist()
        if len(arr) < norg:
            arr += [default] * (norg - len(arr))
        return arr

    cfg["fit"]["energy"] = reset(cfg_energy, default=None)
    cfg["fit"]["energyweight"] = reset(cfg_energyweight)
    cfg["fit"]["energyflag"] = reset(cfg_energyflag)
    cfg["fit"]["energyscatter"] = reset(cfg_energyscatter)

    # Dummy matrix (apparently needed for multi-energy)
    if cfg["attenuators"]["Matrix"][0] == 0 and nenergies > 1:
        cfg["materials"]["Dummy"] = {
            "Comment": "Dummy",
            "CompoundFraction": [1],
            "CompoundList": ["H1"],
            "Density": 1.0,
            "Thickness": 0.0,
        }
        cfg["attenuators"]["Matrix"][0] = 1
        cfg["attenuators"]["Matrix"][1] = "Dummy"
        cfg["attenuators"]["Matrix"][2] = 1.0
        cfg["attenuators"]["Matrix"][3] = 0.0  # thickness in cm


def AdaptPyMcaConfig_mlines(cfg, mlines):

    # Split M-lines
    # /usr/local/lib/python2.7/dist-packages/PyMca5/PyMcaPhysics/xrf/Elements.py
    # /users/opid21/.local/lib/python2.7/site-packages/PyMca5/PyMcaPhysics/xrf/Elements.py
    #
    # You need an adapted pymca version: Elements
    # ElementShellTransitions = [KShell.ElementKShellTransitions,
    #                       KShell.ElementKAlphaTransitions,
    #                       KShell.ElementKBetaTransitions,
    #                       LShell.ElementLShellTransitions,
    #                       LShell.ElementL1ShellTransitions,
    #                       LShell.ElementL2ShellTransitions,
    #                       LShell.ElementL3ShellTransitions,
    #                      [s+"*" for s in MShell.ElementMShellTransitions],
    #                       MShell.ElementM1ShellTransitions,
    #                       MShell.ElementM2ShellTransitions,
    #                       MShell.ElementM3ShellTransitions,
    #                       MShell.ElementM4ShellTransitions,
    #                       MShell.ElementM5ShellTransitions]
    # ElementShellRates = [KShell.ElementKShellRates,
    #                 KShell.ElementKAlphaRates,
    #                 KShell.ElementKBetaRates,
    #                 LShell.ElementLShellRates,
    #                 LShell.ElementL1ShellRates,
    #                 LShell.ElementL2ShellRates,
    #                 LShell.ElementL3ShellRates,
    #                 MShell.ElementMShellRates,
    #                 MShell.ElementM1ShellRates,
    #                 MShell.ElementM2ShellRates,
    #                 MShell.ElementM3ShellRates,
    #                 MShell.ElementM4ShellRates,
    #                 MShell.ElementM5ShellRates]
    # ElementXrays      = ['K xrays', 'Ka xrays', 'Kb xrays', 'L xrays','L1 xrays','L2 xrays','L3 xrays','M xrays','M1 xrays','M2 xrays','M3 xrays','M4 xrays','M5 xrays']
    if "M5 xrays" not in ClassMcaTheory.Elements.ElementXrays:
        msg = "XRF fit: PyMca5.PyMcaPhysics.xrf.Elements is not patched to supported M-line group splitting."
        logger.error(msg)
        raise ImportError(msg)
    for el in mlines:
        if el in cfg["peaks"]:
            if "M" in cfg["peaks"][el]:
                cfg["peaks"][el] = [
                    group for group in cfg["peaks"][el] if group != "M"
                ] + mlines[el]


def AdaptPyMcaConfig_quant(cfg, quant):
    if "flux" in quant:
        cfg["concentrations"]["flux"] = quant["flux"]
    if "time" in quant:
        cfg["concentrations"]["time"] = quant["time"]
    if "area" in quant:
        cfg["concentrations"]["area"] = quant["area"]
    if "distance" in quant:
        cfg["concentrations"]["distance"] = quant["distance"]
    if "anglein" in quant:
        cfg["attenuators"]["Matrix"][4] = quant["anglein"]
    if "angleout" in quant:
        cfg["attenuators"]["Matrix"][5] = quant["angleout"]
    if "anglein" in quant or "angleout" in quant:
        cfg["attenuators"]["Matrix"][7] = (
            cfg["attenuators"]["Matrix"][4] + cfg["attenuators"]["Matrix"][5]
        )


def AdaptPyMcaConfig_fast(cfg):
    if cfg["fit"]["linearfitflag"] == 0:
        cfg["fit"]["linearfitflag"] = 1

    if "strategyflag" not in cfg["fit"]:
        cfg["fit"]["strategyflag"] = 0
    elif cfg["fit"]["strategyflag"]:
        cfg["fit"]["strategyflag"] = 0

    cfg["fit"]["fitweight"] = 0


def AdaptPyMcaConfig_forcebatch(cfg):
    # Force no weights (for spectra with low counts):
    cfg["fit"]["fitweight"] = 0


def AdaptPyMcaConfig_modinfo(cfg, quant, fast):
    ind = instance.asarray(cfg["fit"]["energyflag"]).astype(bool)
    _energy = instance.asarray(cfg["fit"]["energy"])[ind]
    _weights = instance.asarray(cfg["fit"]["energyweight"])[ind]
    _weights = _weights / _weights.sum() * 100
    _scatter = instance.asarray(cfg["fit"]["energyscatter"])[ind]

    info = "\n ".join(
        [
            "{} keV (Rate = {:.2f}%, Scatter {})".format(en, w, "ON" if scat else "OFF")
            for en, w, scat in zip(_energy, _weights, _scatter)
        ]
    )
    if quant:
        info += "\n flux = {:e} s^(-1)\n time = {} s\n active area = {} cm^2\n sample-detector distance = {} cm\n angle IN = {} deg\n angle OUT = {} deg".format(
            cfg["concentrations"]["flux"],
            cfg["concentrations"]["time"],
            cfg["concentrations"]["area"],
            cfg["concentrations"]["distance"],
            cfg["attenuators"]["Matrix"][4],
            cfg["attenuators"]["Matrix"][5],
        )

    if cfg["attenuators"]["Matrix"][0] == 0:
        info += "\n Matrix = None"
    else:
        info += "\n Matrix = {}".format(cfg["attenuators"]["Matrix"][1])
    info += "\n Linear = {}".format("YES" if cfg["fit"]["linearfitflag"] else "NO")
    info += "\n Fast fitting = {}".format("YES" if fast else "NO")
    info += "\n Error propagation = {}".format(
        "Poisson" if cfg["fit"]["fitweight"] else "OFF"
    )
    info += "\n Matrix adjustment = {}".format(
        "ON" if cfg["fit"]["strategyflag"] else "OFF"
    )

    logger.info("XRF fit configuration adapted:\n {}".format(info))


def AdaptPyMcaConfig(cfg, energy, addhigh=0, mlines=None, quant=None, fast=False):
    """
    Args:
        cfg(ConfigDict): pymca configuration
        energy(float): primary beam energy in keV
        addhigh(Optional(num)): add high primary energy with very low weight
        mlines(Optional(dict)): elements (keys) which M line group must be replaced by some M subgroups (values)
        quant(Optional(dict)):
    """
    AdaptPyMcaConfig_energy(cfg, energy, addhigh)
    if mlines:
        AdaptPyMcaConfig_mlines(cfg, mlines)
    if quant and isinstance(quant, dict):
        AdaptPyMcaConfig_quant(cfg, quant)
    if fast:
        AdaptPyMcaConfig_fast(cfg)
    AdaptPyMcaConfig_forcebatch(cfg)
    AdaptPyMcaConfig_modinfo(cfg, quant, fast)


def PerformRoi(filelist, rois, norm=None):
    """ROI XRF spectra in batch with changing primary beam energy.

    Args:
        filelist(list(str)|np.array): spectra to fit
        rois(dict(2-tuple)): ROIs
        norm(Optional(np.array)): normalization array
    Returns:
        dict: {label:nenergies x nfiles,...}
    """
    # Load data
    # Each spectrum (each row) in 1 edf file is acquired at a different energy
    if isinstance(filelist, list):
        dataStack = EDFStack.EDFStack(filelist, dtype=np.float32).data
    else:
        dataStack = filelist
    nfiles, nenergies, nchannels = dataStack.shape

    # Normalization
    if norm is None:
        norm = [1] * nenergies
    else:
        if hasattr(norm, "__iter__"):
            if len(norm) == 1:
                norm = [norm[0]] * nenergies
            elif len(norm) != nenergies:
                raise ValueError(
                    "Expected {} normalization values ({} given)".format(
                        nenergies, len(norm)
                    )
                )
        else:
            norm = [norm] * nenergies

    # ROI
    ret = {}
    for k in rois:
        ret[k] = np.zeros((nenergies, nfiles), dtype=type(dataStack))

    for i in range(nfiles):
        for k, roi in rois.items():
            ret[k][:, i] = np.sum(dataStack[i, :, roi[0] : roi[1]], axis=1) / norm

    return ret


def PerformFit(
    filelist,
    cfgfile,
    energies,
    mlines={},
    norm=None,
    fast=False,
    addhigh=0,
    prog=None,
    plot=False,
):
    """Fit XRF spectra in batch with changing primary beam energy.

    Args:
        filelist(list(str)|np.array): spectra to fit
        cfgfile(str): configuration file to use
        energies(np.array): primary beam energies
        mlines(Optional(dict)): elements (keys) which M line group must be replaced by some M subgroups (values)
        norm(Optional(np.array)): normalization array
        fast(Optional(bool)): fast fitting (linear)
        addhigh(Optional(number)): add higher energy
        prog(Optional(timing.ProgessLogger)): progress object
        plot(Optional(bool))
    Returns:
        dict: {label:nenergies x nfiles,...}
    """

    # Load data
    # Each spectrum (each row) in 1 edf file is acquired at a different energy
    if isinstance(filelist, list):
        dataStack = EDFStack.EDFStack(filelist, dtype=np.float32).data
    else:
        dataStack = filelist

    nfiles, nenergies, nchannels = dataStack.shape

    # MCA channels
    xmin = 0
    xmax = nchannels - 1
    x = np.arange(nchannels, dtype=np.float32)

    # Energies
    if hasattr(energies, "__iter__"):
        if len(energies) == 1:
            energies = [energies[0]] * nenergies
        elif len(energies) != nenergies:
            raise ValueError(
                "Expected {} energies ({} given)".format(nenergies, len(energies))
            )
    else:
        energies = [energies] * nenergies

    # Normalization
    if norm is None:
        norm = [1] * nenergies
    else:
        if hasattr(norm, "__iter__"):
            if len(norm) == 1:
                norm = [norm[0]] * nenergies
            elif len(norm) != nenergies:
                raise ValueError(
                    "Expected {} normalization values ({} given)".format(
                        nenergies, len(norm)
                    )
                )
        else:
            norm = [norm] * nenergies

    # Prepare plot
    if plot:
        fig, ax = plt.subplots()

    # Prepare fit
    # ClassMcaTheory.DEBUG = 1
    mcafit = ClassMcaTheory.McaTheory()
    try:
        mcafit.useFisxEscape(True)
    except:
        pass
    if fast:
        mcafit.enableOptimizedLinearFit()
    else:
        mcafit.disableOptimizedLinearFit()
    cfg = mcafit.configure(ConfigDict.ConfigDict(filelist=cfgfile))

    # Fit at each energy
    if prog is not None:
        prog.setnfine(nenergies * nfiles)

    ret = {}
    for j in range(nenergies):
        # Prepare fit with this energy
        AdaptPyMcaConfig(cfg, energies[j], mlines=mlines, fast=fast, addhigh=addhigh)
        mcafit.configure(cfg)

        # Fit all spectra with this energy
        for i in range(nfiles):
            # Data to fit
            y = dataStack[i, j, :].flatten()
            mcafit.setData(x, y, xmin=xmin, xmax=xmax)

            # Initial parameter estimates
            mcafit.estimate()

            # Fit
            fitresult = mcafit.startfit(digest=0)

            # Extract result
            if plot:
                mcafitresult = mcafit.digestresult()
                ax.cla()

                if (
                    plot == 2
                    or not any(np.isfinite(np.log(mcafitresult["ydata"])))
                    or not any(mcafitresult["ydata"] > 0)
                ):
                    ax.plot(mcafitresult["energy"], mcafitresult["ydata"])
                    ax.plot(mcafitresult["energy"], mcafitresult["yfit"], color="red")
                else:
                    ax.semilogy(mcafitresult["energy"], mcafitresult["ydata"])
                    ax.semilogy(
                        mcafitresult["energy"], mcafitresult["yfit"], color="red"
                    )
                    ax.set_ylim(
                        ymin=np.nanmin(
                            mcafitresult["ydata"][np.nonzero(mcafitresult["ydata"])]
                        )
                    )
                ax.set_title("Primary energy: {} keV".format(energies[j]))
                ax.set_xlabel("Energy (keV)")
                ax.set_ylabel("Intensity (cts)")
                plt.pause(0.0001)
            else:
                mcafitresult = mcafit.imagingDigestResult()

            # Store result
            for k in mcafitresult["groups"]:
                if k not in ret:
                    ret[k] = np.zeros(
                        (nenergies, nfiles), dtype=type(mcafitresult[k]["fitarea"])
                    )
                ret[k][j, i] = mcafitresult[k]["fitarea"] / norm[j]

            if "chisq" not in ret:
                ret["chisq"] = np.zeros((nenergies, nfiles), dtype=type(mcafit.chisq))
            ret["chisq"][j, i] = mcafit.chisq

        # Print progress
        if prog is not None:
            prog.ndonefine(nfiles)
            prog.printprogress()

    return ret


def PerformBatchFit(*args, **kwargs):
    if OutputBuffer is None:
        return PerformBatchFitOld(*args, **kwargs)
    else:
        return PerformBatchFitNew(*args, **kwargs)


def PerformBatchFitHDF5(
    filelist,
    cfg,
    outuri,
    energy=None,
    mlines=None,
    quant=None,
    fast=False,
    addhigh=0,
    **kw
):
    """Fit XRF spectra in batch with one primary beam energy.

        Least-square fitting. If you intend a linear fit, modify the configuration:
          - Get current energy calibration with "Load From Fit"
          - Enable: Perform a Linear Fit
          - Disable: Stripping
          - Strip iterations = 0
        Fast linear least squares:
          - Use SNIP instead of STRIP

    Args:
        filelist(list(str)): spectra to fit
        cfg(str or ConfigDict): configuration file to use
        outuri(h5fs.Path): directory for results
        energy(num): primary beam energy
        mlines(Optional(dict)): elements (keys) which M line group must be replaced by some M subgroups (values)
        fast(Optional(bool)): fast fitting (linear)
        quant(Optional(dict)):
        addhigh(Optional(int)):
    """
    if instance.isstring(cfg):
        cfg = ConfigDict.ConfigDict(filelist=cfg)
    AdaptPyMcaConfig(
        cfg, energy, mlines=mlines, quant=quant, fast=fast, addhigh=addhigh
    )

    # outputDir/outputRoot.h5::/fileEntry/fileProcess
    kw["h5"] = True
    kw["edf"] = False
    kw["outputDir"] = outuri.device.parent.path
    kw["outputRoot"] = os.path.splitext(outuri.device.name)[0]
    kw["fileEntry"] = outuri.parent.path
    kw["fileProcess"] = outuri.name
    outbuffer = OutputBuffer(**kw)
    if fast:
        batch = FastXRFLinearFit.FastXRFLinearFit()
        stack = FastXRFLinearFit.prepareDataStack(filelist)
        kwargs = {
            "y": stack,
            "configuration": cfg,
            "concentrations": bool(quant),
            "refit": 1,
            "outbuffer": outbuffer,
        }
        with outbuffer.saveContext():
            batch.fitMultipleSpectra(**kwargs)
    else:
        split_results = list(zip(*(filename.split("::") for filename in filelist)))
        if len(split_results) == 1:
            selection = None
        else:
            filelist, path_in_file = split_results
            if len(set(path_in_file)) != 1:
                raise ValueError(path_in_file, "HDF5 group must be the same for all")
            filelist = list(filelist)
            selection = {"y": path_in_file[0]}
        kwargs = {
            "filelist": filelist,
            "selection": selection,
            "concentrations": bool(quant),
            "fitfiles": 0,
            "fitconcfile": 0,
            "outbuffer": outbuffer,
        }
        with tempPyMcaConfigFile(cfg) as cfgfilename:
            batch = McaAdvancedFitBatch.McaAdvancedFitBatch(cfgfilename, **kwargs)
            with outbuffer.saveContext():
                batch.processList()


def PerformBatchFitNew(
    filelist,
    outdir,
    outname,
    cfg,
    energy,
    mlines=None,
    quant=None,
    fast=False,
    addhigh=0,
):
    """Fit XRF spectra in batch with one primary beam energy.

        Least-square fitting. If you intend a linear fit, modify the configuration:
          - Get current energy calibration with "Load From Fit"
          - Enable: Perform a Linear Fit
          - Disable: Stripping
          - Strip iterations = 0
        Fast linear least squares:
          - Use SNIP instead of STRIP

    Args:
        filelist(list(str)): spectra to fit
        outdir(str): directory for results
        outname(str): output radix
        cfg(str or ConfigDict): configuration file to use
        energy(num): primary beam energy
        mlines(Optional(dict)): elements (keys) which M line group must be replaced by some M subgroups (values)
        fast(Optional(bool)): fast fitting (linear)
        quant(Optional(dict)):
        addhigh(Optional(int))
    Returns:
        files(list(str)): files produced by pymca
        labels(list(str)): corresponding HDF5 labels
    """
    # Adapt cfg in memory
    if instance.isstring(cfg):
        cfg = ConfigDict.ConfigDict(filelist=cfg)
    AdaptPyMcaConfig(
        cfg, energy, mlines=mlines, quant=quant, fast=fast, addhigh=addhigh
    )
    buncertainties = False
    bconcentrations = bool(quant)

    # Save cfg in temporary file
    outdir = localfs.Path(outdir).mkdir()
    with outdir.temp(name=outname + ".cfg", force=True) as cfgfile:
        cfg.write(cfgfile.path)
        kwargs = {
            "outputDir": outdir.path,
            "fileEntry": outname,
            "h5": False,
            "edf": True,
            "multipage": False,
            "saveFOM": True,
        }
        outbuffer = OutputBuffer(**kwargs)
        if fast:
            batch = FastXRFLinearFit.FastXRFLinearFit()
            stack = FastXRFLinearFit.prepareDataStack(filelist)
            kwargs = {
                "y": stack,
                "configuration": cfg,
                "concentrations": bconcentrations,
                "refit": 1,
                "weight": None,  # None -> from cfg file
                "outbuffer": outbuffer,
            }
        else:
            kwargs = {
                "filelist": filelist,
                "concentrations": bconcentrations,
                "fitfiles": 0,
                "fitconcfile": 0,
                "outbuffer": outbuffer,
            }
            batch = McaAdvancedFitBatch.McaAdvancedFitBatch(cfgfile.path, **kwargs)

        with outbuffer.saveContext():
            if fast:
                batch.fitMultipleSpectra(**kwargs)
            else:
                batch.processList()

    # List of files and labels
    files, labels = [], []
    groups = ["parameters", "massfractions"]
    if buncertainties:
        groups.append("uncertainties")
    for group in groups:
        for label in outbuffer.labels(group, labeltype="filename"):
            filename = outbuffer.filename(".edf", suffix="_" + label)
            labels.append(label)
            files.append(filename)
    if "chisq" in outbuffer:
        labels.append("calc_chisq")
        files.append(outbuffer.filename(".edf", suffix="_chisq"))
    return files, labels


def PerformBatchFitOld(
    filelist,
    outdir,
    outname,
    cfg,
    energy,
    mlines=None,
    quant=None,
    fast=False,
    addhigh=0,
):
    """Fit XRF spectra in batch with one primary beam energy.

        Least-square fitting. If you intend a linear fit, modify the configuration:
          - Get current energy calibration with "Load From Fit"
          - Enable: Perform a Linear Fit
          - Disable: Stripping
          - Strip iterations = 0
        Fast linear least squares:
          - Use SNIP instead of STRIP

    Args:
        filelist(list(str)): spectra to fit
        outdir(str): directory for results
        outname(str): output radix
        cfg(str or ConfigDict): configuration file to use
        energy(num): primary beam energy
        mlines(Optional(dict)): elements (keys) which M line group must be replaced by some M subgroups (values)
        fast(Optional(bool)): fast fitting (linear)
        quant(Optional(dict)):
        addhigh(Optional(int))
    Returns:
        files(list(str)): files produced by pymca
        labels(list(str)): corresponding HDF5 labels
    """
    outdir = localfs.Path(outdir).mkdir()
    if instance.isstring(cfg):
        cfg = ConfigDict.ConfigDict(filelist=cfg)

    with outdir.temp(name=outname + ".cfg", force=True) as cfgfile:
        AdaptPyMcaConfig(
            cfg, energy, mlines=mlines, quant=quant, fast=fast, addhigh=addhigh
        )
        cfg.write(cfgfile.path)

        buncertainties = False
        bconcentrations = bool(quant)
        if fast:
            # Prepare fit
            fastFit = FastXRFLinearFit.FastXRFLinearFit()
            fastFit.setFitConfiguration(cfg)
            dataStack = EDFStack.EDFStack(filelist, dtype=np.float32)

            # Fit
            result = fastFit.fitMultipleSpectra(
                y=dataStack, refit=1, concentrations=bconcentrations
            )

            # Save result and keep filenames + labels
            names = result["names"]
            if bconcentrations:
                names = names[: -len(result["concentrations"])]
            parse = re.compile("^(?P<Z>.+)[_ -](?P<line>.+)$")

            def filename(x):
                return outdir["{}_{}.edf".format(outname, x)].path

            labels = []
            files = []
            j = 0
            for i, name in enumerate(names):
                m = parse.match(name)
                if not m:
                    continue
                m = m.groupdict()
                Z, line = m["Z"], m["line"]

                # Peak area
                label = "{}_{}".format(Z, line)
                f = filename(label)
                edf.saveedf(
                    f, result["parameters"][i], {"Title": label}, overwrite=True
                )
                labels.append(label)
                files.append(f)

                # Error on peak area
                if buncertainties:
                    label = "s{}_{}".format(Z, line)
                    f = filename(label)
                    edf.saveedf(
                        f, result["uncertainties"][i], {"Title": label}, overwrite=True
                    )
                    labels.append(label)
                    files.append(f)

                # Mass fraction
                if bconcentrations and Z.lower() != "scatter":
                    label = "w{}_{}".format(Z, line)
                    f = filename(label)
                    edf.saveedf(
                        f, result["concentrations"][j], {"Title": label}, overwrite=True
                    )
                    labels.append(label)
                    files.append(f)
                    j += 1
        else:
            b = McaAdvancedFitBatch.McaAdvancedFitBatch(
                cfgfile.path,
                filelist=filelist,
                outputdir=outdir.path,
                fitfiles=0,
                concentrations=bconcentrations,
            )
            b.processList()
            filemask = os.path.join(outdir.path, "IMAGES", "*.dat")

            def basename(x):
                return os.path.splitext(os.path.basename(x))[0]

            nbase = len(basename(glob.glob(filemask)[0])) + 1
            filemask = os.path.join(outdir.path, "IMAGES", "*.edf")
            labels = []
            files = []
            for name in sorted(glob.glob(filemask)):
                label = basename(name)[nbase:]
                if label.endswith("mass_fraction"):
                    label = "w" + label[:-14]
                if label == "chisq":
                    label = "calc_chisq"
                labels.append(label)
                files.append(name)
    return files, labels

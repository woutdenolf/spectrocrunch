# -*- coding: utf-8 -*-

import os
import numpy as np
from PyMca5.PyMcaIO import ConfigDict
from ..io import localfs
from ..io import spe
from ..utils import subprocess
from .xrmc import xrmcresult_to_mca


def proc_result(args, out, err, returncode):
    success = returncode == 0
    if not success:
        print("Failed: " + " ".join(args))
        if err:
            print("XMIMSIM errors:")
            print(err)
    return success


def installed():
    return subprocess.installed("xmimsim")


def execute(args, cwd):
    out, err, returncode = subprocess.execute(*args, cwd=cwd, stderr=True)
    if isinstance(out, bytes):
        out = out.decode()
    if isinstance(err, bytes):
        err = err.decode()
    return proc_result(args, out, err, returncode)


def pymcacfg_add_mcinfo(
    configdict,
    outpath,
    p_polarisation=1,
    ninteractions=1,
    multiplicity=100000,
    source_distance=100,
    has_atmosphere=False,
    beamsize=1e-4,
):
    configdict["xrfmc"] = {}
    mcsetup = configdict["xrfmc"]["setup"] = {}
    mcsetup["p_polarisation"] = p_polarisation

    # Point source
    mcsetup["source_diverg_x"] = 0
    mcsetup["source_diverg_y"] = 0
    mcsetup["source_size_x"] = 0
    mcsetup["source_size_y"] = 0
    mcsetup["source_sample_distance"] = source_distance

    # Divergence determined by slits
    # divergence = np.arctan2(beamsize*0.5, source_distance)
    #           = np.arctan2(slit_width*0.5, slit_distance)
    slit_distance = source_distance / 2.0
    slit_width = beamsize / 2.0
    mcsetup["slit_distance"] = slit_distance
    mcsetup["slit_width_x"] = slit_width
    mcsetup["slit_width_y"] = slit_width

    mcsetup["nmax_interaction"] = ninteractions
    # first non-atmospheric layer
    if has_atmosphere:
        mcsetup["layer"] = 2
    else:
        mcsetup["layer"] = 1
    mcsetup["output_dir"] = outpath
    mcsetup["histories"] = multiplicity

    attenuators = configdict["attenuators"]
    for name in "BeamFilter0", "BeamFilter1", "Absorber":
        if name not in attenuators:
            attenuators[name] = [0, "-", 0.0, 0.0, 1.0]


def pymcacfg_to_xmimsimcfg(pymcacfg, xmimsimcfg, **kwargs):
    configdict = ConfigDict.ConfigDict(filelist=pymcacfg)
    outpath = os.path.dirname(xmimsimcfg)
    pymcacfg_add_mcinfo(configdict, outpath, **kwargs)
    configdict.write(xmimsimcfg)


def run_xmimsim_pymca(xmimsimcfg, xmso, pileup=True, escape=True, convolute=True):
    xmsopath = os.path.dirname(xmso)
    xmsofile = os.path.basename(xmso)
    basename = os.path.splitext(xmsofile)[0]
    args = [
        "xmimsim-pymca",
        "--spe-file-unconvoluted={}_lines".format(basename),
        "--verbose",
        "--enable-single-run",
    ]
    if pileup:
        pileup = "--enable-pile-up"
    else:
        pileup = "--disable-pile-up"
    args.append(pileup)
    if escape:
        escape = "--enable-escape-peaks"
    else:
        escape = "--disable-escape-peaks"
    args.append(escape)
    if convolute:
        args.append("--spe-file={}_convoluted".format(basename))
    args += [xmimsimcfg, xmsofile]
    return execute(args, xmsopath)


def xmso_to_xmsi(xmso, xmsi):
    xmsipath = os.path.dirname(xmsi)
    args = ["xmso2xmsi", xmso, xmsi]
    if execute(args, xmsipath):
        patch_xmsi(xmso, xmsi)
        return True
    else:
        return False


def patch_xmsi(xmso, xmsi):
    with open(xmsi) as f:
        content = f.readlines()
    content = [x.rstrip() for x in content]
    try:
        i = content.index(r"    <outputfile/>")
        content[i] = r"    <outputfile>{}</outputfile>".format(xmso)
    except ValueError:
        pass
    try:
        i = content.index(r"    <pulse_width>0</pulse_width>")
        content[i] = r"    <pulse_width>1e-12</pulse_width>"
    except ValueError:
        pass
    try:
        i = [r"1e-5</pulse_width>" in line for line in content].index(True)
        content[i] = content[i].replace(r"1e-5</pulse_width>", r"1e-12</pulse_width>")
    except ValueError:
        pass
    with open(xmsi, mode="w") as f:
        f.write("\n".join(content))


def xmsi_to_xrmc(xmsi, xrmcpath, pileup=True):
    args = ["xmsi2xrmc", xmsi]
    if pileup:
        pileup = "--enable-pile-up"
    else:
        pileup = "--disable-pile-up"
    args.append(pileup)
    return execute(args, xrmcpath)


def run(
    outpath,
    pymcahandle=None,
    pymcacfg=None,
    pileup=True,
    escape=True,
    convolute=True,
    outradix="out",
    runxrmc=False,
    **kwargs
):
    """
    Args:
        outpath(str):
        pymcahandle(PyMcaHandle): overwrites all other arguments
        pymcacfg(str)
        pileup(bool)
        escape(bool)
        outradix(str)
        runxrmc(bool): for comparison
        convolute(bool)
        **kwargs: see `pymcacfg_add_mcinfo`

    Returns:
        success(bool)
    """
    outpath = localfs.Path(outpath).mkdir()
    if pymcahandle is not None:
        pymcacfg = str(outpath["{}.cfg".format(outradix)])
        pymcahandle.savepymca(pymcacfg)
        pileup = pymcahandle.pileup
        escape = pymcahandle.escape
        kwargs["ninteractions"] = pymcahandle.ninteractions
    else:
        pymcacfg = str(localfs.Path(pymcacfg).copy(outpath[pymcacfg]))
    xmimsimcfg = str(outpath["{}_xmimsim.cfg".format(outradix)])
    pymcacfg_to_xmimsimcfg(pymcacfg, xmimsimcfg, **kwargs)
    xmso = str(outpath["{}.xmso".format(outradix)])
    if not run_xmimsim_pymca(
        xmimsimcfg, xmso, escape=escape, pileup=pileup, convolute=convolute
    ):
        return False
    xmsi = str(outpath["{}.xmsi".format(outradix)])
    if not xmso_to_xmsi(xmso, xmsi):
        return False
    if runxrmc:
        xrmcpath = outpath["xrmc"].mkdir()
        if not xmsi_to_xrmc(xmsi, str(xrmcpath), pileup=pileup):
            return False
        if not execute(["xrmc", "input.dat"], str(xrmcpath)):
            return False
        xrmcfile = str(xrmcpath["convoluted_spectra.dat"])
        mcafile = str(xrmcpath["{}.mca".format(outradix)])
        xrmcresult_to_mca(xrmcfile, mcafile, mode="w")
    return True


def loadxmimsimresult(outpath, outradix="out", convoluted=False):
    outpath = localfs.Path(outpath)
    if convoluted:
        suffix = "convoluted"
    else:
        suffix = "lines"
    fmt = "{}_{}_{{}}.spe".format(outradix, suffix)
    i = 1
    while outpath[fmt.format(i)].exists:
        i += 1
    i -= 1
    mca, channels, energy, coeff = spe.read(str(outpath[fmt.format(i)]))
    zero, gain = coeff
    info = {"xenergy": energy, "zero": zero, "gain": gain}
    return mca, info

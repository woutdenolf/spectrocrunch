# -*- coding: utf-8 -*-

import os
from glob import glob
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from PyMca5.PyMca import ArraySave
from PyMca5.PyMcaCore import XiaEdf

from ..io.spec import spec
from ..xrf.parse_xia import parse_xia_esrf
from ..xrf.fit import PerformFit as fitter
from ..xrf.fit import PerformRoi as roisummer
from ..utils.timing import ProgressLogger

import pylab
import logging

# def angletoenergy(angle,dspacing):
#    """
#    Args:
#        angle: degrees
#        dspacing: monochromator d-spacing
#    Returns:
#        energy (keV)
#    """
#    hc = 4.13566743E-8 * 299792458
#    return hc/(2*dspacing*np.sin(angle*np.pi/180))


def processNotSynchronized(
    specfile,
    specnumbers,
    destpath,
    detectorcfg,
    mlines={},
    replacebasedir=None,
    showelement=None,
    dtcor=True,
    fastfitting=True,
    energyshift=0,
    plot=False,
    bkgxas=0,
    bkgflux=0,
    normmedian=False,
    rois=None,
    counters=None,
    energylabel="arr_energyM",
    iodetlabel="arr_iodet",
    timelabel="arr_mtime",
):
    """
    XRF fitting of XANES spectra (fit repeats separately and add interpolated results because no energy synchronization)

    Args:
        specfile(str): name of the spec file
        specnumbers(list(list(int))): list of lists of spec numbers
        destpath(str): directory for saving the result
        detectorcfg(list(str)): config files for fitting (one per detector)
        replacebasedir(Optional(2-tuple)): replace first with second in the data directory extracted from the spec file
        dtcor(Optional(True)): correct spectrum for deadtime before fitting
        fastfitting(Optional(True)): linear fitting or non-linear
        showelement(Optional(str)): element to be plotted
        energyshift(Optional(num)): energy shift in keV
        plot(Optional(bool)): plot results
        mlines(Optional(dict)): elements (keys) which M line group must be replaced by some M subgroups (values)
        bkgxas(Optional(Num)): subtract from XAS spectrum
        bkgflux(Optional(Num)): subtract from iodet signal (cts/sec)
        normmedian(Optional(bool)): normalize the XRF count normalization
        rois(Optional(list(dict(2-tuple)))): ROIs instead of fitting
        counters(Optional(dict)): list of counters to be treated as the XRF counts
    """
    addbeforefit = (
        True  # refers to multiple detectors, repeats are always added afterwards
    )

    # Open spec file
    sf = spec(specfile)

    # Prepare
    nxasspectra = len(specnumbers)
    nrepeats = [len(l) for l in specnumbers]
    nxasspectraT = sum(nrepeats)
    if dtcor:
        parsename = "dtcor"
    else:
        parsename = "copy"
    logger = logging.getLogger(__name__)
    prog = ProgressLogger(logger)
    if not hasattr(detectorcfg, "__iter__"):
        detectorcfg = [detectorcfg]
    if isinstance(rois, dict):
        rois = [rois]
    if counters is None:
        counters = {}
    ncounters = len(counters)
    counteroutnames = counters.keys()
    counterinnames = [counters[c]["name"] for c in counters]
    counterbkg = [counters[c]["bkg"] for c in counters]
    counterminlog = [counters[c]["minlog"] for c in counters]

    # Loop over spectra
    off = 0
    prog.start()
    for i in range(nxasspectra):
        # XAS spectrum: sum of all repeats and detectors
        xasspectrum = {}

        # Loop over repeats
        nrepeats = len(specnumbers[i])
        for j in range(nrepeats):
            # Get energy and iodet
            has_timelabel = bool(timelabel)
            if not has_timelabel:
                timelabel = iodetlabel
            data, info = sf.getdata(
                specnumbers[i][j], [energylabel, iodetlabel, timelabel] + counterinnames
            )
            realtime = sf.scancommand(specnumbers[i][j])["time"]
            if not has_timelabel:
                data[:, 2] = realtime

            data[:, 0] += energyshift
            energyj = data[:, 0][:, np.newaxis]
            norm = (data[:, 1] / realtime - bkgflux) * data[:, 2]

            if ncounters > 0:
                ctrs = data[:, 2:]
                for c in range(ncounters):
                    ctrs[:, c] = (ctrs[:, c] / realtime - counterbkg[c]) * data[:, 2]

            if normmedian:
                norm /= np.median(norm)

            # Parse xia files
            datadir = info["DIRECTORY"]
            if len(replacebasedir) == 2:
                datadir = datadir.replace(replacebasedir[0], replacebasedir[1])
                detectorcfg = [
                    f.replace(replacebasedir[0], replacebasedir[1]) for f in detectorcfg
                ]

            scanname = info["RADIX"]
            scannumber = int(info["ZAP SCAN NUMBER"])
            if dtcor:
                parsename = "dtcor"
            else:
                parsename = "copy"
            parsename = (
                "%%0%dd_%s"
                % (np.int(np.floor(np.log10(nxasspectraT))) + 1, parsename)
                % (off)
            )
            if j == 0:
                destradix = scanname
                outdir = os.path.join(destpath, destradix + "_data")
            filestofit, detnums = parse_xia_esrf(
                datadir,
                scanname,
                scannumber,
                outdir,
                parsename,
                deadtime=dtcor,
                add=addbeforefit,
            )
            ndets = len(filestofit)

            # Intialize progress counter
            if i == 0 and j == 0:
                prog.setn(nxasspectraT * ndets)

            # Fit, normalize and add spectra from all detector
            xasspectrumj = {}
            for k in range(ndets):
                idet = detnums[k]

                if len(filestofit[k]) != 0:

                    if rois is None:
                        if len(detectorcfg) == 1:
                            cfg = detectorcfg[0]
                        else:
                            cfg = detectorcfg[k]

                        # Perform fitting
                        xasresults = fitter(
                            filestofit[k],
                            cfg,
                            energyj,
                            mlines=mlines,
                            norm=norm,
                            fast=fastfitting,
                            prog=prog,
                            plot=plot,
                        )
                    else:
                        if len(rois) == 1:
                            roisk = rois[0]
                        else:
                            roisk = rois[k]

                        # Perform ROI summing
                        xasresults = roisummer(filestofit[k], roisk, norm=norm)

                    if len(xasspectrumj) == 0:
                        xasspectrumj = xasresults
                    # elif energy_ref is None:
                    else:
                        for group in xasresults:
                            xasspectrumj[group] += xasresults[group]

            # Add normalized sum of counters
            for c in range(ncounters):
                tmp = ctrs[:, c]
                if counterminlog[c]:
                    tmp = -np.log(tmp)
                xasspectrumj[counteroutnames[c]] = tmp[:, np.newaxis]

            # Add this repeat to the previous repeats (if any)
            if len(xasspectrum) == 0:
                xasspectrum = xasspectrumj
                energy = energyj
            else:
                for group in xasspectrumj:
                    spl = InterpolatedUnivariateSpline(
                        energyj, xasspectrumj[group], ext=0
                    )
                    xasspectrum[group] += spl(energy)

            # Show
            if showelement in xasspectrum and plot:
                pylab.clf()
                pylab.plot(energy, xasspectrum[showelement][:, 0])
                pylab.title(
                    "Spec #{}: {}/I0 (Summed repeats = {})".format(
                        specnumbers[i][j], showelement, j + 1
                    )
                )
                pylab.pause(0.01)

            # Show progress
            prog.ndone(ndets)
            prog.printprogress()

        # What we want:
        #    XAS1 = sum(I)/sum(I0) = mu(fluo).rho.d
        # What we have:
        #    XAS1 = sum(I/I0) = nrepeats.mu(fluo).rho.d
        for group in xasspectrum:
            xasspectrum[group] /= nrepeats

        # Save XAS spectrum (for each element)
        if specnumbers[i][0] == specnumbers[i][-1]:
            outname = "{}_{:03d}".format(destradix, specnumbers[i][0])
        else:
            outname = "{}_{:03d}_{:03d}.sum".format(
                destradix, specnumbers[i][0], specnumbers[i][-1]
            )
        fileName = os.path.join(destpath, outname + ".dat")
        if not os.path.exists(destpath):
            os.makedirs(destpath)

        xasspectrum["energy"] = energy
        labels = [k.replace(" ", "-") for k in xasspectrum]
        ArraySave.save2DArrayListAsASCII(
            list(xasspectrum.values()), fileName, labels=labels
        )
        logger.info("Saved XAS spectrum {}.".format(fileName))


def processEnergySynchronized(
    specfile,
    specnumbers,
    destpath,
    pymcacfg,
    mlines={},
    replacebasedir=(),
    showelement=None,
    dtcor=True,
    fastfitting=True,
    energyshift=0,
    plot=False,
    bkgxas=0,
    bkgflux=0,
    normmedian=False,
    rois=None,
    counters=None,
    energylabel="arr_energyM",
    iodetlabel="arr_iodet",
    timelabel="arr_mtime",
):
    """
    XRF fitting of XANES spectra (add spectra from repeats because of energy synchronization)

    Args:
        specfile(str): name of the spec file
        specnumbers(list(list(int))): list of lists of spec numbers
        destpath(str): directory for saving the result
        pymcacfg(str): config file for fitting
        replacebasedir(Optional(2-tuple)): replace first with second in the data directory extracted from the spec file
        dtcor(Optional(True)): correct spectrum for deadtime before fitting
        fastfitting(Optional(True)): linear fitting or non-linear
        showelement(Optional(str)): element to be plotted
        energyshift(Optional(num)): energy shift in keV
        plot(Optional(bool)): plot results
        mlines(Optional(dict)): elements (keys) which M line group must be replaced by some M subgroups (values)
        bkgxas(Optional(Num)): subtract from XAS spectrum
        bkgflux(Optional(Num)): subtract from iodet signal (cts/sec)
        normmedian(Optional(bool)): normalize the XRF count normalization
        rois(Optional(dict(2-tuple))): ROIs instead of fitting
        counters(Optional(dict)): list of counters to be treated as the XRF counts
    """

    # /users/blissadm/spec/macros/zap/zapxmap.mac
    #  xmap_x1_00 gets its counts from (see ZAP_PSEUDO): _zap_xmap_roi_calc
    #  xmap_x1_00 = sum(XRF[a:b])/mtime*time
    #
    # /users/blissadm/spec/macros/zap/zaptools.mac
    #  arr_iodet = iodet/arr_mtime*realtime
    #
    # -> xanes = sum(XRF[a:b])/(mtime*arr_iodet)*time

    # Open spec file
    sf = spec(specfile)

    # Prepare
    nxasspectra = len(specnumbers)
    nrepeats = [len(l) for l in specnumbers]
    nxasspectraT = sum(nrepeats)
    logger = logging.getLogger(__name__)
    prog = ProgressLogger(logger)
    if not os.path.exists(destpath):
        os.makedirs(destpath)
    if counters is None:
        counters = {}
    ncounters = len(counters)
    counteroutnames = counters.keys()
    counterinnames = [counters[c]["name"] for c in counters]
    counterbkg = [counters[c]["bkg"] for c in counters]

    # Loop over spectra
    prog.setn(nxasspectra)
    prog.start()
    for i in range(nxasspectra):
        # XAS spectrum: sum of all repeats and detectors
        xasspectrum = {}
        nrepeats = len(specnumbers[i])
        xrfinfo = [{}] * nrepeats

        # Get spec info
        for j in range(nrepeats):
            # Get energy and iodet
            has_timelabel = bool(timelabel)
            if not has_timelabel:
                timelabel = iodetlabel
            data, info = sf.getdata(
                specnumbers[i][j], [energylabel, iodetlabel, timelabel] + counterinnames
            )
            realtime = sf.scancommand(specnumbers[i][j])["time"]
            if not has_timelabel:
                data[:, 2] = realtime

            data[:, 0] += energyshift

            # Check energy synchronization
            if "energy" in xasspectrum:
                if xasspectrum["nenergy"] != data.shape[0]:
                    raise ValueError(
                        "Number of energies in spec scan {} are not the same as for spec scan {}".format(
                            specnumbers[i], specnumbers[0]
                        )
                    )
                if not np.allclose(
                    xasspectrum["energy"], data[:, 0], rtol=0, atol=1e-6
                ):
                    raise ValueError(
                        "Energies in spec scan {} are not synchronized with energies in {}".format(
                            specnumbers[i], specnumbers[0]
                        )
                    )
            else:
                xasspectrum["nenergy"] = data.shape[0]
                xasspectrum["energy"] = data[:, 0]
                xasspectrum["norm"] = np.empty(
                    (data.shape[0], nrepeats), dtype=data.dtype
                )
                if ncounters != 0:
                    for c in counterinnames:
                        xasspectrum["counter_" + c] = np.empty(
                            (data.shape[0], nrepeats), dtype=data.dtype
                        )

            # data[:,1] = iodet/data[:,2]*realtime
            # norm = (data[:,1]/realtime - bkgflux)*data[:,2]
            #      = iodet - bkgflux*data[:,2]
            xasspectrum["norm"][:, j] = (data[:, 1] / realtime - bkgflux) * data[:, 2]
            xrfinfo[j] = info
            if ncounters != 0:
                for c in range(ncounters):
                    xasspectrum["counter_" + counterinnames[c]][:, j] = (
                        data[:, 3 + c] / realtime - counterbkg[c]
                    ) * data[:, 2]

        if normmedian:
            xasspectrum["norm"] /= np.median(xasspectrum["norm"])

        # Generate normalized XRF spectra to be fitted
        for j in range(nrepeats):

            norm = xasspectrum["norm"][:, j][:, np.newaxis]
            info = xrfinfo[j]

            # Parse xia files
            datadir = info["DIRECTORY"]
            if len(replacebasedir) == 2:
                datadir = datadir.replace(replacebasedir[0], replacebasedir[1])
                if pymcacfg is not None:
                    pymcacfg = pymcacfg.replace(replacebasedir[0], replacebasedir[1])
            scanname = info["RADIX"]
            scannumber = int(info["ZAP SCAN NUMBER"])

            # Sum, dt correction and I0 normalize
            fs = os.path.join(
                datadir, "%s_xia[0-9]*_%04d_0000_*.edf" % (scanname, scannumber)
            )
            detfiles = sorted(glob(fs))
            if len(detfiles) == 0:
                logger.warning("No files found with filter {}".format(fs))
            fs = os.path.join(
                datadir, "%s_xiast_%04d_0000_*.edf" % (scanname, scannumber)
            )
            stfile = glob(fs)
            if len(stfile) == 0:
                logger.warning("No files found with filter {}".format(fs))

            if len(detfiles) == 0:
                xia = None
                if "data" not in xasspectrum:
                    xasspectrum["data"] = norm * 0
            else:
                xia = XiaEdf.XiaEdfScanFile(stfile[0], detfiles)
                err = xia.sum(deadtime=dtcor)
                if "data" in xasspectrum:
                    xasspectrum["data"] += xia.data / norm
                else:
                    xasspectrum["data"] = xia.data / norm

        # What we want:
        #    XAS = sum(I)/sum(I0) = mu(fluo).rho.d
        # What we have:
        #    XAS = sum(I/I0) = nrepeats.mu(fluo).rho.d
        xasspectrum["data"] /= nrepeats

        # Subtract background
        xasspectrum["data"] -= bkgxas
        xasspectrum["data"][xasspectrum["data"] < 0] = 0

        # Fit spectrum or take ROI
        energy = xasspectrum["energy"]
        if specnumbers[i][0] == specnumbers[i][-1]:
            outname = "{}_{:03d}".format(scanname, specnumbers[i][0])
        else:
            outname = "{}_{:03d}_{:03d}.sum".format(
                scanname, specnumbers[i][0], specnumbers[i][-1]
            )

        if rois is None and xia is not None:
            # Save XRF spectra to be fitted (not needed, just for checking the fit afterwards)
            fileName = os.path.join(destpath, outname + ".edf")
            xia.data = xasspectrum["data"]
            xia.save(fileName, 1)

            # Fit xas spectrum
            datastack = xasspectrum["data"][np.newaxis, ...]
            xasresults = fitter(
                datastack,
                pymcacfg,
                energy,
                mlines=mlines,
                fast=fastfitting,
                prog=prog,
                plot=plot,
            )

            # Show fit result
            if showelement in xasresults and plot:
                pylab.clf()
                pylab.plot(energy, xasresults[showelement][:, 0])
                pylab.title(
                    "Spec #{}-#{}: {}/I0 ({} repeats)".format(
                        specnumbers[i][0], specnumbers[i][-1], showelement, nrepeats
                    )
                )
                pylab.pause(0.01)
        else:
            # datastack = xasspectrum["data"][np.newaxis,...]
            # xasresults = roisummer(datastack,rois)

            # More straightforward:
            xasresults = {}
            if rois is not None:
                nen, nchan = xasspectrum["data"].shape
                for label, roi in rois.items():
                    xasresults[label] = np.sum(
                        xasspectrum["data"][:, roi[0] : roi[1]], axis=1
                    )[:, None]

        # if True:
        #    import matplotlib.pyplot as plt
        #    plt.plot(xasresults["Fe-Ka"][:,0],label="script")
        #    plt.plot(data[:,-1],label="arr_absorp3")
        #    plt.title("{}: #{}".format(specfile,specnumbers[0][0]))
        #    plt.legend()
        #    plt.show()
        #    exit()

        # Add energy to result
        xasresults["energy"] = energy[:, np.newaxis]

        # Add normalized sum of counters
        for c in counters:
            tmp = np.sum(
                xasspectrum["counter_" + counters[c]["name"]] / xasspectrum["norm"],
                axis=1,
            )
            tmp /= nrepeats

            if counters[c]["minlog"]:
                tmp = -np.log(tmp)

            xasresults[c] = tmp[:, np.newaxis]

        # Add norm
        xasresults["norm"] = xasspectrum["norm"]

        # Save XAS spectrum (for each element)
        fileName = os.path.join(destpath, outname + ".dat")
        labels = [k.replace(" ", "-") for k in xasresults]
        ArraySave.save2DArrayListAsASCII(
            list(xasresults.values()), fileName, labels=labels
        )
        logger.info("Saved XAS spectrum {}.".format(fileName))

        # Show progress
        prog.ndone(1)
        prog.printprogress()

from PyMca5.PyMcaCore import XiaCorrect
from glob import glob
import os


def log_dummy(message, verbose_level=None, verbose_ask=None):
    pass


def parse_xia_esrf(
    datadir,
    scanname,
    scannumber,
    outdir,
    outname,
    exclude_detectors=[],
    deadtime=True,
    add=True,
    willbefitted=True,
    deadtimeifsingle=True,
):
    """Parse an XIA XRF scan with ESRF naming convention
        Input files:
            [scanname]_xia[xianum]_[scannumber]_0000_[linenumber].edf
            xianum = detector number or "st"
        Output files:
            [scanname]_[outname]_xia[xianum]_[scannumber]_0000_[linenumber].edf
            xianum = detector number or "S1" for sum

    Args:
        datadir(list(str)): directory of the xia files
        scanname(list(str)): radix of the xia files
        scannumber(list(int)): scan number of the xia files
        outdir(str): directory for the corrected/summed xia files if any
        outname(str): radix for the corrected/summed xia files if any
        exclude_detectors(Optional(list[int])): detector numbers to be excluded
        deadtime(Optional(bool)): apply deadtime correction to spectra
        add(Optional(bool)): add spectra (after deadtime correction if any)
        willbefitted(Optional(bool)): spectra will be fitted
        deadtimeifsingle(Optional(bool)): deadtime can be skipped for a single detector (will be done afterwards)
    """

    # Raw data
    if not isinstance(datadir, list):
        datadir = [datadir]
    if not isinstance(scanname, list):
        scanname = [scanname]
    if not isinstance(scannumber, list):
        scannumber = [scannumber]
    if not isinstance(scannumber, list):
        scannumber = [scannumber]

    files_raw = []
    for ddir, sname, snum in zip(datadir, scanname, scannumber):
        filemask = os.path.join(ddir, "%s_xia*_%04d_0000_*.edf" % (sname, snum))
        files_raw.extend(glob(filemask))
    files_raw = sorted(files_raw)

    # Check files
    xiafiles = XiaCorrect.parseFiles(files_raw, log_cb=log_dummy)
    # - files which only differ by their xianum are grouped together (i.e. each linescan or ct)
    # - missing but existing xianum's are added
    if xiafiles is None:
        return [[]], []

    # Exclude detectors
    if len(exclude_detectors) != 0:
        xiafiles = [
            [f for f in fs if f.getDetector() not in exclude_detectors]
            for fs in xiafiles
        ]

    # Detector numbers
    tmp = list(exclude_detectors)
    tmp.append(None)  # exclude "st" detector
    detnums = [f.getDetector() for f in xiafiles[0] if f.getDetector() not in tmp]
    ndet = len(detnums)

    # Only count detectors and not save sum/dtcorrection?
    if ndet < 2:
        add = False
        if not deadtimeifsingle:
            deadtime = False

    if ndet == 0:
        add = False
        deadtime = False
        willbefitted = False

    if not willbefitted and not deadtime and not add:
        return [[]], detnums

    # DT correction and optional summation
    if add or deadtime:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        if add:
            sums = [detnums]
        else:
            sums = None

        XiaCorrect.correctFiles(
            xiafiles,
            deadtime=deadtime,
            livetime=0,
            sums=sums,
            avgflag=0,
            outdir=outdir,
            outname=outname,
            force=1,
            log_cb=None,
        )

        # Output files
        files_out = []
        for ddir, sname, snum in zip(datadir, scanname, scannumber):
            if add:
                filemask = os.path.join(
                    outdir, "%s_%s_xiaS1_%04d_0000_*.edf" % (ddir, sname, snum)
                )
            else:
                filemask = os.path.join(
                    outdir, "%s_%s_xia[0-9]*_%04d_0000_*.edf" % (ddir, sname, snum)
                )
            files_out.extend(glob(filemask))
        files_out = sorted(files_out)
    else:
        files_out = sorted(
            [
                file.get()
                for detfiles in xiafiles
                for file in detfiles
                if not file.isStat()
            ]
        )

    # Group according to detector
    if add:
        files_out = [files_out]
    else:
        nf = len(files_out)
        nfdet = nf // ndet
        files_out = [files_out[i : i + nfdet] for i in range(0, nf, nfdet)]

    return files_out, detnums

# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from PyMca5.PyMcaCore import XiaCorrect
import glob
import os

def log_dummy(message, verbose_level=None, verbose_ask=None):
    pass

def parse_xia_esrf(datadir,scanname,scannumber,outdir,outname,exclude_detectors=[],deadtime=True,add=True):
    """Parse an XIA XRF scan with ESRF naming
        Input files:
            [scanname]_xia[xianum]_[scannumber]_0000_[linenumber].edf
            xianum = detector number or "st"
        Output files:
            [scanname]_[outname]_xia[xianum]_[scannumber]_0000_[linenumber].edf
            xianum = detector number or "S1" for sum
    """

    # Raw data
    filemask = os.path.join(datadir,"%s_xia*_%04d_0000_*.edf"%(scanname,scannumber))
    files_raw = sorted(glob.glob(filemask))
    
    # Check files
    xiafiles = XiaCorrect.parseFiles(files_raw) # list of list
    if xiafiles is None:
        return [[]],[]

    # Exclude detectors
    if len(exclude_detectors)!=0:
        xiafiles = [[f for f in fs if f.getDetector() not in exclude_detectors] for fs in xiafiles]

    # Detector numbers
    tmp = list(exclude_detectors)
    tmp.append(None) # exclude "st" detector
    detnums = [f.getDetector() for f in xiafiles[0] if f.getDetector() not in tmp]
    ndet = len(detnums)
    if ndet <= 0:
        return [[]],detnums
    if ndet == 1:
        add = False

    # DT correction and optional summation
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    sums = [detnums] if add else None
    XiaCorrect.correctFiles(xiafiles, deadtime=deadtime, livetime=0, sums=sums, 
                            avgflag=0, outdir=outdir, outname=outname, force=1,
                            log_cb=None)

    # Output files
    if add:
        filemask = os.path.join(outdir,"%s_%s_xiaS1_%04d_0000_*.edf"%(scanname,outname,scannumber))
    else:
        filemask = os.path.join(outdir,"%s_%s_xia[0-9]*_%04d_0000_*.edf"%(scanname,outname,scannumber))
    files_out = sorted(glob.glob(filemask))

    # Group according to detector
    if add:
        files_out = [files_out]
    else:
        nf = len(files_out)
        nfdet = nf//ndet
        files_out = [files_out[i:i+nfdet] for i in range(0,nf,nfdet)]

    return files_out,detnums


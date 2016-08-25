# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
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

import os
from glob import glob
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from PyMca5.PyMca import ArraySave
from PyMca5.PyMcaCore import XiaEdf

from spectrocrunch.io.spec import spec
from spectrocrunch.xrf.parse_xia import parse_xia_esrf
from spectrocrunch.xrf.fit import PerformFit as fitter
from spectrocrunch.common.timing import progress

import pylab
import logging

#def angletoenergy(angle,dspacing):
#    """
#    Args:
#        angle: degrees
#        dspacing: monochromator d-spacing
#    Returns:
#        energy (keV)
#    """
#    hc = 4.13566743E-8 * 299792458
#    return hc/(2*dspacing*np.sin(angle*np.pi/180))

def processNotSynchronized(specfile,specnumbers,destpath,detectorcfg,mlines={},replacebasedir=None,showelement=None,dtcor=True,fastfitting=True,energyshift=0,plot=False):
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
    """

    energylabel = 'arr_energyM'
    iodetlabel = 'arr_iodet'
    addbeforefit = True # refers to multiple detectors, repeats are always added afterwards
  
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
    prog = progress(logger)
    if not hasattr(detectorcfg,"__iter__"):
        detectorcfg = [detectorcfg]

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
            data,info = sf.getdata(specnumbers[i][j],[energylabel,iodetlabel])
            data[:,0] += energyshift
            energyj = data[:,0].reshape(data.shape[0],1)
            norm = data[:,1]

            # Parse xia files
            datadir = info["DIRECTORY"]
            if len(replacebasedir) == 2:
                datadir = datadir.replace(replacebasedir[0],replacebasedir[1])
                
            scanname = info["RADIX"]
            scannumber = int(info["ZAP SCAN NUMBER"])
            if dtcor:
                parsename = "dtcor"
            else:
                parsename = "copy"
            parsename = "%%0%dd_%s"%(np.int(np.floor(np.log10(nxasspectraT)))+1,parsename)%(off)
            if j==0:
                destradix = scanname
                outdir = os.path.join(destpath,destradix+"_data")
            filestofit,detnums = parse_xia_esrf(datadir,scanname,scannumber,outdir,parsename,deadtime=dtcor,add=addbeforefit)
            ndets = len(filestofit)

            # Intialize progress counter
            if i==0 and j==0:
                prog.setn(nxasspectraT*ndets)

            # Fit, normalize and add spectra from all detector
            xasspectrumj = {}
            for k in range(ndets):
                idet = detnums[k]

                if len(filestofit[k])!= 0:
                    if len(detectorcfg)==1:
                        cfg = detectorcfg[0]
                    else:
                        cfg = detectorcfg[k]

                    # Perform fitting
                    fitresults = fitter(filestofit[k],cfg,energyj,mlines=mlines,norm=norm,fast=fastfitting,prog=prog,plot=plot)
                    
                    if len(xasspectrumj)==0:
                        xasspectrumj = fitresults
                    elif energy_ref is None:
                        for group in fitresults:
                            xasspectrumj[group] += fitresults[group]

            # Add this repeat to the previous repeats (if any)
            if len(xasspectrum)==0:
                xasspectrum = xasspectrumj
                energy = energyj
            else:
                for group in xasspectrumj:
                    spl = InterpolatedUnivariateSpline(energyj,xasspectrumj[group],ext=0)
                    xasspectrum[group] += spl(energy)

            # Show 
            if showelement in xasspectrum and plot:
                pylab.clf()
                pylab.plot(energy,xasspectrum[showelement][:,0])
                pylab.title("Spec #{}: {}/I0 (Summed repeats = {})".format(specnumbers[i][j],showelement,j+1))
                pylab.pause(0.01)

            # Show progress
            prog.ndone(ndets)
            prog.printprogress()

        # Save XAS spectrum (for each element)
        outname = destradix+'_'+'_'.join([str(d) for d in specnumbers[i]])
        fileName = os.path.join(destpath, outname+".dat")
        if not os.path.exists(destpath):
            os.makedirs(destpath)

        xasspectrum["energy"] = energy
        labels = [k.replace(' ','-') for k in xasspectrum]
        ArraySave.save2DArrayListAsASCII(xasspectrum.values(), fileName, labels=labels)
        logger.info("Saved XAS spectrum {}.".format(fileName))

def processEnergySynchronized(specfile,specnumbers,destpath,pymcacfg,mlines={},replacebasedir=None,showelement=None,dtcor=True,fastfitting=True,energyshift=0,plot=False):
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
    """

    energylabel = 'arr_energyM'
    iodetlabel = 'arr_iodet'
  
    # Open spec file
    sf = spec(specfile)

    # Prepare
    nxasspectra = len(specnumbers)
    nrepeats = [len(l) for l in specnumbers]
    nxasspectraT = sum(nrepeats)
    logger = logging.getLogger(__name__)
    prog = progress(logger)
    if not os.path.exists(destpath):
        os.makedirs(destpath)

    # Loop over spectra
    prog.setn(nxasspectra)
    prog.start()
    for i in range(nxasspectra):
        # XAS spectrum: sum of all repeats and detectors
        xasspectrum = {}
        nrepeats = len(specnumbers[i])

        # Get spec info
        for j in range(nrepeats):
            # Get energy and iodet
            data,info = sf.getdata(specnumbers[i][j],[energylabel,iodetlabel])
            data[:,0] += energyshift

            # Check energy synchronization
            if "energy" in xasspectrum:
                if xasspectrum["nenergy"]!=data.shape[0]:
                    raise ValueError("Number of energies in spec scan {} is are not the same as for spec scan {}".format(specnumbers[i],specnumbers[0]))
                if not np.allclose(xasspectrum["energy"],data[:,0],rtol=0,atol=1e-6):
                    raise ValueError("Energies in spec scan {} is are not synchronized with energies in {}".format(specnumbers[i],specnumbers[0]))
            else:
                xasspectrum["nenergy"] = data.shape[0]
                xasspectrum["energy"] = data[:,0]
                xasspectrum["norm"] = np.empty((data.shape[0],nrepeats),dtype=data.dtype)
            xasspectrum["norm"][:,j] = data[:,1]
        
        xasspectrum["norm"] /= np.median(xasspectrum["norm"])
 
        # Generate XRF spectra to be fitted
        for j in range(nrepeats):
            norm = xasspectrum["norm"][:,j].reshape((xasspectrum["nenergy"],1))

            # Parse xia files
            datadir = info["DIRECTORY"]
            if len(replacebasedir) == 2:
                datadir = datadir.replace(replacebasedir[0],replacebasedir[1])
            scanname = info["RADIX"]
            scannumber = int(info["ZAP SCAN NUMBER"])

            # Sum, dt correction and I0 normalize
            fs = os.path.join(datadir,"%s_xia[0-9]*_%04d_0000_*.edf"%(scanname,scannumber))
            detfiles = sorted(glob(fs))
            if len(detfiles)==0:
                logger.error("No files found with filter {}".format(fs))
            fs = os.path.join(datadir,"%s_xiast_%04d_0000_*.edf"%(scanname,scannumber))
            stfile = glob(fs)
            if len(stfile)==0:
                logger.error("No files found with filter {}".format(fs))
            xia = XiaEdf.XiaEdfScanFile(stfile[0], detfiles)
            err = xia.sum(deadtime=dtcor)
            if "data" in xasspectrum:
                xasspectrum["data"] += xia.data/norm
            else:
                xasspectrum["data"] = xia.data/norm

        # Save XRF spectra to be fitted
        outname = scanname+'_'+'_'.join([str(d) for d in specnumbers[i]])
        fileName = os.path.join(destpath, outname+".edf")
        xia.data = xasspectrum["data"]
        xia.save(fileName, 1)

        # Fit xanes spectrum
        datastack = xasspectrum["data"].reshape((1,)+xasspectrum["data"].shape)
        energy = xasspectrum["energy"]
        fitresults = fitter(datastack,pymcacfg,energy,mlines=mlines,fast=fastfitting,prog=prog,plot=plot)
        fitresults["energyM"] = energy.reshape(xasspectrum["nenergy"],1)
 
        # Show
        if showelement in fitresults and plot:
            pylab.clf()
            pylab.plot(energy,fitresults[showelement][:,0])
            pylab.title("Spec #{}-#{}: {}/I0 ({} repeats)".format(specnumbers[i][0],specnumbers[i][-1],showelement,nrepeats))
            pylab.pause(0.01)

        # Save XAS spectrum (for each element)
        fileName = os.path.join(destpath, outname+".dat")
        labels = [k.replace(' ','-') for k in fitresults]
        ArraySave.save2DArrayListAsASCII(fitresults.values(), fileName, labels=labels)
        logger.info("Saved XAS spectrum {}.".format(fileName))

        # Show progress
        prog.ndone(1)
        prog.printprogress()


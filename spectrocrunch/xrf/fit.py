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

from PyMca5.PyMcaPhysics.xrf import McaAdvancedFitBatch
from PyMca5.PyMcaPhysics.xrf import FastXRFLinearFit as FastXRFLinearFitBase
from PyMca5.PyMcaPhysics.xrf import ClassMcaTheory
from PyMca5.PyMca import EDFStack
from PyMca5.PyMca import ArraySave
from PyMca5.PyMcaIO import EdfFile
from PyMca5.PyMcaIO import ConfigDict

import numpy as np
import re
import os
import glob
import collections
import matplotlib.pyplot as plt
#import warnings
#warnings.filterwarnings("ignore")
import logging

from ..common import instance
from ..io import edf

logger = logging.getLogger(__name__)


class FastXRFLinearFit(FastXRFLinearFitBase.FastXRFLinearFit):

    def fitMultipleSpectra(self,concentrations=False,**kwargs):
        result = super(FastXRFLinearFit,self).fitMultipleSpectra(concentrations=concentrations,**kwargs)
        return result
        
        config = self._mcaTheory.getConfiguration()
    

def ReadPyMcaConfigFile(filename):
    # Read the configuration
    if not os.path.exists(filename):
        raise IOError("File <%s> does not exists" % filename)
    cfg = ConfigDict.ConfigDict()
    cfg.read(filename)
    if len(cfg)==0:
        raise IOError("File <%s> couldn't be loaded" % filename)
    return cfg

def AdaptPyMcaConfigFile(filename,energy,addhigh=True,mlines={},quant={}):
    cfg = ReadPyMcaConfigFile(filename)
    AdaptPyMcaConfig(cfg,energy,addhigh=addhigh,mlines=mlines,quant=quant)
    cfg.write(filename)

def AdaptPyMcaConfig(cfg,energy,addhigh=True,mlines={},quant={}):
    """
    Args:
        cfg(ConfigDict): pymca configuration
        energy(float): primary beam energy in keV
        addhigh(Optional(num)): add high primary energy with very low weight
        mlines(Optional(dict)): elements (keys) which M line group must be replaced by some M subgroups (values)
        quant(Optional(dict)): 
    """

    # Nothing to do
    if np.isnan(energy) and not mlines and not quant:
        return

    if not np.isnan(energy):
        # Adapt the cfg
        
        tmp = cfg["fit"]["energyweight"]
        if instance.isarray(tmp):
            ftype = type(tmp[0])
        else:
            ftype = type(tmp)

        tmp = cfg["fit"]["energyflag"]
        if instance.isarray(tmp):
            itype = type(tmp[0])
        else:
            itype = type(tmp)
        
        n = 1+addhigh

        # Adapt energy
        sourcelines = [None]*n
        sourcelines[0] = ftype(energy)
        if addhigh:
            sourcelines[1] = ftype(10*energy)
        cfg["fit"]["energy"] = sourcelines

        sourcelines = [ftype(0)]*n
        if addhigh:
            sourcelines[0] = ftype(1e100)
            sourcelines[1] = ftype(1e-5)
        else:
            sourcelines[0] = ftype(1)
        cfg["fit"]["energyweight"] = sourcelines

        sourcelines = [itype(0)]*n
        sourcelines[0] = itype(1)
        if addhigh:
            sourcelines[1] = itype(1)
        cfg["fit"]["energyflag"] = sourcelines

        sourcelines = [itype(0)]*n
        sourcelines[0] = itype(1)
        cfg["fit"]["energyscatter"] = sourcelines

        # Dummy matrix (aparently needed for multi-energy)
        if (cfg["attenuators"]["Matrix"][0]==0 and addhigh):
            density = cfg["materials"]["Air"]["Density"]
            cfg["attenuators"]["Matrix"][0] = 1
            cfg["attenuators"]["Matrix"][1] = "Air"
            cfg["attenuators"]["Matrix"][2] = density
            cfg["attenuators"]["Matrix"][3] = density*0 # thickness in cm

    # Split M-lines
    # /usr/local/lib/python2.7/dist-packages/PyMca5/PyMcaPhysics/xrf/Elements.py
    #
    # You need an adapted pymca version: Elements
    #ElementShellTransitions = [KShell.ElementKShellTransitions,
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
    #ElementShellRates = [KShell.ElementKShellRates,
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
    #ElementXrays      = ['K xrays', 'Ka xrays', 'Kb xrays', 'L xrays','L1 xrays','L2 xrays','L3 xrays','M xrays','M1 xrays','M2 xrays','M3 xrays','M4 xrays','M5 xrays']

    if mlines:
        if "M5 xrays" not in ClassMcaTheory.Elements.ElementXrays:
            msg = "XRF fit: PyMca5.PyMcaPhysics.xrf.Elements is not patched to supported M-line group splitting."
            logger.error(msg)
            raise ImportError(msg)
        for el in mlines:
            if el in cfg["peaks"]:
                if "M" in cfg["peaks"][el]:
                    cfg["peaks"][el] = [group for group in cfg["peaks"][el] if group != "M"] + mlines[el]

    if quant:
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
            cfg["attenuators"]["Matrix"][7] = cfg["attenuators"]["Matrix"][4]+cfg["attenuators"]["Matrix"][5]

    # Show info
    _energy = instance.asarray(cfg["fit"]["energy"])
    _weights = instance.asarray(cfg["fit"]["energyweight"])
    _weights = _weights/_weights.sum()*100
    sourceinfo = "\n ".join(["{} keV({:.2f}%)".format(en,w) for en,w in zip(_energy,_weights)])
    if quant:
        fluxinfo = "\n flux = {:e} s^(-1)\n time = {} s\n active area = {} cm^2\n sample-detector distance = {} cm\n angle IN = {} deg\n angle OUT = {} deg".\
                format(cfg["concentrations"]["flux"],\
                       cfg["concentrations"]["time"],\
                       cfg["concentrations"]["area"],\
                       cfg["concentrations"]["distance"],\
                       cfg["attenuators"]["Matrix"][4],\
                       cfg["attenuators"]["Matrix"][5])
    else:
        fluxinfo = "\n"
        
    logger.info("XRF fit configuration adapted:\n {}{}".format(sourceinfo,fluxinfo))


def PerformRoi(filelist,rois,norm=None):
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
    if isinstance(filelist,list):
        dataStack = EDFStack.EDFStack(filelist, dtype=np.float32).data
    else:
        dataStack = filelist
    nfiles,nenergies,nchannels = dataStack.shape

    # Normalization
    if norm is None:
        norm = [1]*nenergies
    else:
        if hasattr(norm,"__iter__"):
            if len(norm)==1:
                norm = [norm[0]]*nenergies
            elif len(norm)!=nenergies:
                raise ValueError("Expected {} normalization values ({} given)".format(nenergies,len(norm)))
        else:
            norm = [norm]*nenergies

    # ROI
    ret = {}
    for k in rois:
        ret[k] = np.zeros((nenergies,nfiles),dtype=type(dataStack))

    for i in range(nfiles):
        for k,roi in rois.items():
            ret[k][:,i] = np.sum(dataStack[i,:,roi[0]:roi[1]],axis=1)/norm

    return ret

def PerformFit(filelist,cfgfile,energies,mlines={},norm=None,fast=True,prog=None,plot=False):
    """Fit XRF spectra in batch with changing primary beam energy.

    Args:
        filelist(list(str)|np.array): spectra to fit
        cfgfile(str): configuration file to use
        energies(np.array): primary beam energies
        mlines(Optional(dict)): elements (keys) which M line group must be replaced by some M subgroups (values)
        norm(Optional(np.array)): normalization array
        fast(Optional(bool)): fast fitting (linear)
        prog(Optional(timing.progress)): progress object
        plot(Optional(bool))
    Returns:
        dict: {label:nenergies x nfiles,...}
    """

    # Load data
    # Each spectrum (each row) in 1 edf file is acquired at a different energy
    if isinstance(filelist,list):
        dataStack = EDFStack.EDFStack(filelist, dtype=np.float32).data
    else:
        dataStack = filelist

    nfiles,nenergies,nchannels = dataStack.shape

    # MCA channels
    xmin = 0
    xmax = nchannels-1
    x = np.arange(nchannels, dtype=np.float32)

    # Energies
    if hasattr(energies,"__iter__"):
        if len(energies)==1:
            energies = [energies[0]]*nenergies
        elif len(energies)!=nenergies:
            raise ValueError("Expected {} energies ({} given)".format(nenergies,len(energies)))
    else:
        energies = [energies]*nenergies

    # Normalization
    if norm is None:
        norm = [1]*nenergies
    else:
        if hasattr(norm,"__iter__"):
            if len(norm)==1:
                norm = [norm[0]]*nenergies
            elif len(norm)!=nenergies:
                raise ValueError("Expected {} normalization values ({} given)".format(nenergies,len(norm)))
        else:
            norm = [norm]*nenergies

    # Prepare plot
    if plot:
        fig, ax = plt.subplots()

    # Prepare fit
    #ClassMcaTheory.DEBUG = 1
    mcafit = ClassMcaTheory.McaTheory()
    try:
        mcafit.useFisxEscape(True)
    except:
        pass
    if fast:
        mcafit.enableOptimizedLinearFit()
    else:
        mcafit.disableOptimizedLinearFit()
    cfg = mcafit.configure(ReadPyMcaConfigFile(cfgfile))

    # Fit at each energy
    if prog is not None:
        prog.setnfine(nenergies*nfiles)

    ret = {}
    for j in range(nenergies):
        # Prepare fit with this energy
        AdaptPyMcaConfig(cfg,energies[j],mlines=mlines)
        mcafit.configure(cfg)

        # Fit all spectra with this energy
        for i in range(nfiles):
            # Data to fit
            y = dataStack[i,j,:].flatten()
            mcafit.setData(x,y,xmin=xmin,xmax=xmax)

            # Initial parameter estimates
            mcafit.estimate()

            # Fit
            fitresult = mcafit.startfit(digest=0)

            # Extract result
            if plot:
                mcafitresult = mcafit.digestresult()
                ax.cla()

                if plot==2 or not any(np.isfinite(np.log(mcafitresult["ydata"]))) or not any(mcafitresult["ydata"]>0):
                    ax.plot(mcafitresult["energy"],mcafitresult["ydata"])
                    ax.plot(mcafitresult["energy"],mcafitresult["yfit"],color='red')
                else:
                    ax.semilogy(mcafitresult["energy"],mcafitresult["ydata"])
                    ax.semilogy(mcafitresult["energy"],mcafitresult["yfit"],color='red')
                    ax.set_ylim(ymin=np.nanmin(mcafitresult["ydata"][np.nonzero(mcafitresult["ydata"])]))
                ax.set_title("Primary energy: {} keV".format(energies[j]))
                ax.set_xlabel("Energy (keV)")
                ax.set_ylabel("Intensity (cts)")
                plt.pause(0.0001)
            else:
                mcafitresult = mcafit.imagingDigestResult()

            # Store result
            for k in mcafitresult["groups"]:
                if k not in ret:
                    ret[k] = np.zeros((nenergies,nfiles),dtype=type(mcafitresult[k]["fitarea"]))
                ret[k][j,i] = mcafitresult[k]["fitarea"]/norm[j]

            if "chisq" not in ret:
                ret["chisq"] = np.zeros((nenergies,nfiles),dtype=type(mcafit.chisq))
            ret["chisq"][j,i] = mcafit.chisq

        # Print progress
        if prog is not None:
            prog.ndonefine(nfiles)
            prog.printprogress()

    return ret

def PerformBatchFit(filelist,outdir,outname,cfgfile,energy,mlines={},quant={},fast=True):
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
        cfgfile(str): configuration file to use
        energy(num): primary beam energy
        mlines(Optional(dict)): elements (keys) which M line group must be replaced by some M subgroups (values)
        fast(Optional(bool)): fast fitting (linear)
        quant(Optional(dict)): 
    """

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Adapt file (not adapting the fitobject's member variables because it's unclear 
    # what other things need to be changed when changing the energy)
    AdaptPyMcaConfigFile(cfgfile,energy,mlines=mlines,quant=quant)

    if fast:
        # Prepare fit
        fastFit = FastXRFLinearFit()
        fastFit.setFitConfigurationFile(cfgfile)
        dataStack = EDFStack.EDFStack(filelist, dtype=np.float32)

        # Fit
        buncertainties = False
        bconcentrations = bool(quant)
        result = fastFit.fitMultipleSpectra(y=dataStack,weight=0,refit=1,concentrations=bconcentrations)

        # Save result and keep filenames + labels
        names = result['names']
        if bconcentrations:
            names = names[:-len(result["concentrations"])]

        parse = re.compile("^(?P<Z>.+)[_ -](?P<line>.+)$")
        filename = lambda x: os.path.join(outdir,"{}_{}.edf".format(outname,x))
        out = collections.OrderedDict()
        for i,name in enumerate(names):
            m = parse.match(name)
            if not m:
                continue

            m = parse.match(name).groupdict()
            Z,line = m["Z"],m["line"]
            
            # Peak area
            label = "{}-{}".format(Z,line)
            f = filename("{}_{}".format(Z,line))
            edf.saveedf(f,\
                        result['parameters'][i],\
                        {'Title': label},overwrite=True)
            out[label] = f
            
            # Error on peak area
            if buncertainties:
                label = "s({}-{})".format(Z,line)
                f = filename("s({}_{})".format(Z,line))
                edf.saveedf(f,\
                            result['uncertainties'][i],\
                            {'Title': label},overwrite=True)
                out[label] = f
                
            # Mass fraction
            if bconcentrations and Z.lower()!="scatter":
                label = "w({}-{})".format(Z,line)
                f = filename("C({}_{})".format(Z,line))
                edf.saveedf(f,\
                            result['concentrations'][i],\
                            {'Title': label},overwrite=True)
                out[label] = f
                
        labels = out.keys()
        files = out.values()
    else:
        # Parallelize this:
        b = McaAdvancedFitBatch.McaAdvancedFitBatch(cfgfile,filelist=filelist,outputdir=outdir,fitfiles=0)
        b.processList()

        #TODO: process filenames like in fast mode
        filemask = os.path.join(outdir,"IMAGES","*.edf")
        files = sorted(glob.glob(filemask))
        files = [f for f in files if "chisq" not in f]
        labels = ["_".join(os.path.splitext(os.path.basename(f))[0].split("_")[-2:]) for f in files]
    
    return files,labels





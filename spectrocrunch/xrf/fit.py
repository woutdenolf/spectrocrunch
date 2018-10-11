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
from PyMca5.PyMcaPhysics.xrf import FastXRFLinearFit
from PyMca5.PyMcaPhysics.xrf import ClassMcaTheory
from PyMca5.PyMca import EDFStack
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

from ..utils import instance
from ..io import edf
from ..io import utils as ioutils

logger = logging.getLogger(__name__)

def ReadPyMcaConfigFile(filename):
    # Read the configuration
    if not os.path.exists(filename):
        raise IOError("File <%s> does not exists" % filename)
    cfg = ConfigDict.ConfigDict()
    cfg.read(filename)
    if len(cfg)==0:
        raise IOError("File <%s> couldn't be loaded" % filename)
    return cfg

def AdaptPyMcaConfigFile(filename,*args,**kwargs):
    cfg = ReadPyMcaConfigFile(filename)
    AdaptPyMcaConfig(cfg,*args,**kwargs)
    cfg.write(filename)

def AdaptPyMcaConfig_energy(cfg,energy,addhigh):
    if not np.isfinite(energy):
        return

    ind = instance.asarray(cfg["fit"]["energyflag"]).astype(bool)
    norg = len(ind)
    nenergies = ind.sum()+addhigh
            
    def extract(name,default=np.nan):
        arr = cfg["fit"][name]
        if instance.isarray(arr):
            arr = [instance.asnumber(v) for v in arr]
        arr = instance.asarray(arr)
        
        # Select based on energyflag
        narr = len(arr)
        if narr<norg:
            arr = np.append(arr,[default]*(norg-narr))
        arr = arr[0:norg][ind]
        
        # At least nenergies
        narr = len(arr)
        if narr<nenergies:
            arr = np.append(arr,[default]*(nenergies-narr))
        return arr
        
    cfg_energy = extract("energy",default=np.nan)
    cfg_energyweight = extract("energyweight",default=np.nan)
    cfg_energyflag = extract("energyflag",default=1)
    cfg_energyscatter = extract("energyscatter",default=0)
    
    cfg_energy = cfg_energy/cfg_energy[0]*energy
    cfg_energyweight = cfg_energyweight/cfg_energyweight[0]
    
    for i in range(nenergies):
        if not np.isfinite(cfg_energy[i]):
            if i==0:
                cfg_energy[i] = energy
            else:
                cfg_energy[i] = 10*energy
        if not np.isfinite(cfg_energyweight[i]):
            if i==0:
                cfg_energyweight[i] = 1
            else:
                cfg_energyweight[i] = 1e-10

    def reset(arr,default=0):
        arr = arr.tolist()
        if len(arr)<norg:
            arr += [default]*(norg-len(arr))
        return arr

    cfg["fit"]["energy"] = reset(cfg_energy,default=None)
    cfg["fit"]["energyweight"] = reset(cfg_energyweight)
    cfg["fit"]["energyflag"] = reset(cfg_energyflag)
    cfg["fit"]["energyscatter"] = reset(cfg_energyscatter)
    
    # Dummy matrix (apparently needed for multi-energy)
    if (cfg["attenuators"]["Matrix"][0]==0 and nenergies>1):
        cfg["materials"]["Dummy"] = {'Comment': 'Dummy', 'CompoundFraction': [1], 'CompoundList': ['H1'], 'Density': 1.0, 'Thickness': 0.0}
        cfg["attenuators"]["Matrix"][0] = 1
        cfg["attenuators"]["Matrix"][1] = "Dummy"
        cfg["attenuators"]["Matrix"][2] = 1.0
        cfg["attenuators"]["Matrix"][3] = 0. # thickness in cm

def AdaptPyMcaConfig_mlines(cfg):
        
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

    if "M5 xrays" not in ClassMcaTheory.Elements.ElementXrays:
        msg = "XRF fit: PyMca5.PyMcaPhysics.xrf.Elements is not patched to supported M-line group splitting."
        logger.error(msg)
        raise ImportError(msg)
    for el in mlines:
        if el in cfg["peaks"]:
            if "M" in cfg["peaks"][el]:
                cfg["peaks"][el] = [group for group in cfg["peaks"][el] if group != "M"] + mlines[el]

def AdaptPyMcaConfig_quant(cfg,quant):
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

def AdaptPyMcaConfig_fast(cfg):
    if cfg["fit"]["linearfitflag"]==0:
        cfg["fit"]["linearfitflag"] = 1
    
    if "strategyflag" not in cfg["fit"]:
        cfg["fit"]["strategyflag"] = 0
    elif cfg["fit"]["strategyflag"]:
        cfg["fit"]["strategyflag"] = 0

    cfg['fit']['fitweight'] = 0 # Bug in pymca?

def AdaptPyMcaConfig_forcebatch(cfg):
    # Force no weights (for spectra with low counts): 
    cfg['fit']['fitweight'] = 0

def AdaptPyMcaConfig_modinfo(cfg,quant):
    ind = instance.asarray(cfg["fit"]["energyflag"]).astype(bool)
    _energy = instance.asarray(cfg["fit"]["energy"])[ind]
    _weights = instance.asarray(cfg["fit"]["energyweight"])[ind]
    _weights = _weights/_weights.sum()*100
    _scatter = instance.asarray(cfg["fit"]["energyscatter"])[ind]
    
    info = "\n ".join(["{} keV (Rate = {:.2f}%, Scatter {})".format(en,w,"ON" if scat else "OFF") for en,w,scat in zip(_energy,_weights,_scatter)])
    if quant:
        info += "\n flux = {:e} s^(-1)\n time = {} s\n active area = {} cm^2\n sample-detector distance = {} cm\n angle IN = {} deg\n angle OUT = {} deg".\
                format(cfg["concentrations"]["flux"],\
                       cfg["concentrations"]["time"],\
                       cfg["concentrations"]["area"],\
                       cfg["concentrations"]["distance"],\
                       cfg["attenuators"]["Matrix"][4],\
                       cfg["attenuators"]["Matrix"][5])
    
    if cfg["attenuators"]["Matrix"][0]==0:
        info += "\n Matrix = None"
    else:
        info += "\n Matrix = {}".format(cfg["attenuators"]["Matrix"][1])
    info += "\n Linear = {}".format("YES" if cfg["fit"]["linearfitflag"] else "NO")
    info += "\n Error propagation = {}".format("Poisson" if cfg['fit']['fitweight'] else "OFF")
    info += "\n Matrix adjustment = {}".format("ON" if cfg["fit"]["strategyflag"] else "OFF")
    
    logger.info("XRF fit configuration adapted:\n {}".format(info))
    
def AdaptPyMcaConfig(cfg,energy,addhigh=True,mlines=None,quant=None,fast=False):
    """
    Args:
        cfg(ConfigDict): pymca configuration
        energy(float): primary beam energy in keV
        addhigh(Optional(num)): add high primary energy with very low weight
        mlines(Optional(dict)): elements (keys) which M line group must be replaced by some M subgroups (values)
        quant(Optional(dict)): 
    """
    AdaptPyMcaConfig_energy(cfg,energy,addhigh)
    if mlines:
        AdaptPyMcaConfig_mlines(cfg)
    if quant:
        AdaptPyMcaConfig_quant(cfg,quant)
    if fast:
        AdaptPyMcaConfig_fast(cfg)
    AdaptPyMcaConfig_forcebatch(cfg)
    
    AdaptPyMcaConfig_modinfo(cfg,quant)
    
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

def PerformFit(filelist,cfgfile,energies,mlines={},norm=None,fast=False,prog=None,plot=False):
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
        AdaptPyMcaConfig(cfg,energies[j],mlines=mlines,fast=fast)
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

def PerformBatchFit(filelist,outdir,outname,cfgfile,energy,mlines=None,quant=None,fast=False):
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
    Returns:
        files(list(str)): files produced by pymca
        labels(list(str)): corresponding HDF5 labels
    """

    # Adapt file
    ioutils.mkdir(outdir)

    with ioutils.Copy(cfgfile,os.path.join(outdir,outname+".cfg")) as cfgfile:
        AdaptPyMcaConfigFile(cfgfile,energy,mlines=mlines,quant=quant,fast=fast)

        buncertainties = False
        bconcentrations = bool(quant)
            
        if fast:
            # Prepare fit
            fastFit = FastXRFLinearFit.FastXRFLinearFit()
            fastFit.setFitConfigurationFile(cfgfile)
            dataStack = EDFStack.EDFStack(filelist, dtype=np.float32)

            # Fit
            result = fastFit.fitMultipleSpectra(y=dataStack,refit=1,concentrations=bconcentrations)

            # Save result and keep filenames + labels
            names = result['names']
            if bconcentrations:
                names = names[:-len(result["concentrations"])]

            parse = re.compile("^(?P<Z>.+)[_ -](?P<line>.+)$")
            filename = lambda x: os.path.join(outdir,"{}_{}.edf".format(outname,x))
            labels = []
            files = []
            
            j = 0
            for i,name in enumerate(names):
                m = parse.match(name)
                if not m:
                    continue
                m = m.groupdict()
                Z,line = m["Z"],m["line"]
                
                # Peak area
                label = "{}-{}".format(Z,line)
                f = filename("{}_{}".format(Z,line))
                edf.saveedf(f,\
                            result['parameters'][i],\
                            {'Title': label},overwrite=True)
                labels.append(label)
                files.append(f)
                
                # Error on peak area
                if buncertainties:
                    label = "s{}-{}".format(Z,line)
                    f = filename("s{}_{}".format(Z,line))
                    edf.saveedf(f,\
                                result['uncertainties'][i],\
                                {'Title': label},overwrite=True)
                    labels.append(label)
                    files.append(f)
                    
                # Mass fraction
                if bconcentrations and Z.lower()!="scatter":
                    label = "w{}-{}".format(Z,line)
                    f = filename("w{}_{}".format(Z,line))
                    edf.saveedf(f,\
                                result['concentrations'][j],\
                                {'Title': label},overwrite=True)
                    labels.append(label)
                    files.append(f)
                    j += 1
        else:
            # TODO: parallelize this
            b = McaAdvancedFitBatch.McaAdvancedFitBatch(cfgfile,filelist=filelist,outputdir=outdir,
                                                        fitfiles=0,concentrations=bconcentrations)
            b.processList()

            filemask = os.path.join(outdir,"IMAGES","*.dat")
            basename = lambda x: os.path.splitext(os.path.basename(x))[0]
            nbase = len(basename(glob.glob(filemask)[0]))+1

            filemask = os.path.join(outdir,"IMAGES","*.edf")
            labels = []
            files = []
            
            for name in sorted(glob.glob(filemask)):
                label = basename(name)[nbase:]
                if label.endswith("mass_fraction"):
                    label = "w"+label[:-14]
                label = label.replace("_","-")
                if label=="chisq":
                    label = "calc_chisq"
                labels.append(label)
                files.append(name)

        return files,labels


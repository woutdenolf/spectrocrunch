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

from PyMca5.PyMcaCore import XiaEdf
from PyMca5.PyMcaPhysics.xrf import McaAdvancedFitBatch
from PyMca5.PyMcaPhysics.xrf import FastXRFLinearFit
from PyMca5.PyMca import EDFStack
from PyMca5.PyMca import ArraySave
from PyMca5.PyMcaIO import EdfFile
from PyMca5.PyMcaIO import ConfigDict

import numpy as np
import re
import os
import glob

import warnings
warnings.filterwarnings("ignore")

def AdaptEnergy(cfg,energy):
    # Read the configuration
    if not os.path.exists(cfg):
        raise IOError("File <%s> does not exists" % cfg)
    configuration = ConfigDict.ConfigDict()
    configuration.read(cfg)

    # Adapt the configuration
    ftype = type(configuration["fit"]["energyweight"][0])
    itype = type(configuration["fit"]["energyflag"][0])
    n = len(configuration["fit"]["energy"])

    # Adapt energy
    sourcelines = [None]*n
    sourcelines[0] = ftype(energy)
    sourcelines[1] = ftype(2*energy)
    configuration["fit"]["energy"] = sourcelines

    sourcelines = [ftype(0)]*n
    sourcelines[0] = ftype(1-1e-100)
    sourcelines[1] = ftype(1e-100)
    configuration["fit"]["energyweight"] = sourcelines

    sourcelines = [itype(0)]*n
    sourcelines[0] = itype(1)
    sourcelines[1] = itype(1)
    configuration["fit"]["energyflag"] = sourcelines

    sourcelines = [itype(0)]*n
    sourcelines[0] = itype(1)
    configuration["fit"]["energyscatter"] = sourcelines

    # Dummy matrix (aparently needed for multi-energy)
    if (configuration["attenuators"]["Matrix"][0]==0):
        density = configuration["materials"]["Air"]["Density"]
        configuration["attenuators"]["Matrix"][0] = 1
        configuration["attenuators"]["Matrix"][1] = "Air"
        configuration["attenuators"]["Matrix"][2] = density
        configuration["attenuators"]["Matrix"][3] = density*0

    # Write the configuration
    configuration.write(cfg)

def PerformBatchFit(filelist,outdir,outname,cfg,energy,fast=True):
    # Least-square fitting. If you intend a linear fit, modify the configuration:
    #   - Get current energy calibration with "Load From Fit"
    #   - Enable: Perform a Linear Fit
    #   - Disable: Stripping
    #   - Strip iterations = 0
    # Fast linear least squares:
    #   - Use SNIP instead of STRIP

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Adapt file (not adapting the fitobject's member variables because it's unclear 
    # what other things need to be changed when changing the energy)
    if energy is not np.nan:
        AdaptEnergy(cfg,energy)

    if fast:
        # Prepare fit
        fastFit = FastXRFLinearFit.FastXRFLinearFit()
        fastFit.setFitConfigurationFile(cfg)
        dataStack = EDFStack.EDFStack(filelist, dtype=np.float32)

        # Fit
        result = fastFit.fitMultipleSpectra(y=dataStack,weight=0,refit=1,concentrations=0)
        buncertainties = False

        # Results to save
        images = result['parameters']
        imageNames = result['names']

        nImages = images.shape[0]
        n = nImages
        if buncertainties:
            n += len(result['uncertainties'])
        imageList = [None] * n
        labels = [None] * n
        j = 0
        for i in range(nImages):
            name = imageNames[i].replace(" ","-")
            labels[j] = name
            imageList[j] = images[i]
            j += 1
            if not imageNames[i].startswith("C(") and buncertainties:
                # fitted parameter
                labels[j] = "s(%s)" % name
                imageList[j] = result['uncertainties'][i]
                j += 1

        # Save result (similar to save2DArrayListAsEDF)
        #prefix = outname #+ "_" + re.split("[_.]+",filelist[0])[-2] + "_to_" + re.split("[_.]+",filelist[-1])[-2]
        n = len(imageList)
        files = [None]*n
        for i in range(n):
            filename = os.path.join(outdir,outname+"_"+labels[i].replace("-","_")+".edf")
            if os.path.exists(filename):
                try:
                    os.remove(filename)
                except OSError:
                    pass
            edfout = EdfFile.EdfFile(filename, access="ab")
            edfout.WriteImage({'Title': labels[i]},
                                    imageList[i], Append=1)
            del edfout  # force file close
            files[i] = filename

        #fileName = os.path.join(outdir, outname+".dat")
        #ArraySave.save2DArrayListAsASCII(imageList, fileName, labels=labels)
    else:
        # Parallelize this:
        b = McaAdvancedFitBatch.McaAdvancedFitBatch(cfg,filelist=filelist,outputdir=outdir,fitfiles=0)
        b.processList()

        filemask = os.path.join(outdir,"IMAGES","*.edf")
        files = sorted(glob.glob(filemask))
        files = [f for f in files if "chisq" not in f]
        labels = ["_".join(os.path.splitext(os.path.basename(f))[0].split("_")[-2:]) for f in files]
    
    return files,labels


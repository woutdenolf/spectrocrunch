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
import os
import numpy as np
import h5py
import logging

from . import batch
from .fluoxas import process
from ..xrf import create_hdf5_imagestacks
from ..utils import instance

logger = logging.getLogger(__name__)

def fluoxas(samplename,datasetname,scannumbers,mapnumbers,cfgfiles,**parameters):
    if len(scannumbers)!=len(mapnumbers):
        raise RuntimeError("fluoXAS map numbers must be equal to the mapnumbers")
    jobname = batch.jobname("multi",(samplename,datasetname,scannumbers,mapnumbers,cfgfiles),parameters)
    
    instrument = create_hdf5_imagestacks.getinstrument(parameters)
    radix,subdir = instrument.xrflocation(samplename,datasetname,type="dynamic")

    sourcepath = [os.path.join(parameters["proposaldir"],subdir,"{}_fluoXAS_{}".format(radix,nr)) for nr in scannumbers]
    scanname = ["{}_fluoXAS_{}".format(radix,nr) for nr in scannumbers]
    out =  "{}_fmap{}".format(radix,scannumbers[0])

    processdata(jobname,sourcepath,scanname,mapnumbers,cfgfiles,out,fluoxas=True,**parameters)

def multi(samplename,datasetname,mapnumbers,cfgfiles,**parameters):
    jobname = batch.jobname("multi",(samplename,datasetname,mapnumbers,cfgfiles),parameters)
    
    instrument = create_hdf5_imagestacks.getinstrument(parameters)
    radix,subdir = instrument.xrflocation(samplename,datasetname,type="dynamic")
    sourcepath = [os.path.join(parameters["proposaldir"],subdir)]

    scanname = [radix]
    scannumbers = [mapnumbers]
    out =  "{}_map{}_{}".format(radix,mapnumbers[0],mapnumbers[-1])

    processdata(jobname,sourcepath,scanname,scannumbers,cfgfiles,out,multi=True,**parameters)

def single(samplename,datasetname,mapnumber,cfgfiles,**parameters):
    jobname = batch.jobname("single",(samplename,datasetname,mapnumber,cfgfiles),parameters)
    
    instrument = create_hdf5_imagestacks.getinstrument(parameters)
    radix,subdir = instrument.xrflocation(samplename,datasetname,type="dynamic")
    sourcepath = [os.path.join(parameters["proposaldir"],subdir)]

    scanname = [radix]
    scannumbers = [[mapnumber]]
    out =  "{}_map{}".format(radix,mapnumber)

    processdata(jobname,sourcepath,scanname,scannumbers,cfgfiles,out,**parameters)
    
def processdata(jobname,*args,**kwargs):
    if "jobs" in kwargs:
        kwargs["jobs"].append((jobname,processdata_exec,args,kwargs))
    else:
        processdata_exec(*args,**kwargs)
        
def processdata_exec(sourcepath,scanname,scannumbers,cfgfiles,out,
                    fluoxas=False,multi=False,geometry=None,resultsdir=None,
                    resultssubdir='crunched',refimageindex=None,**kwargs):
    parameters = dict(kwargs)

    # Basic input
    parameters['sourcepath'] = sourcepath
    parameters['scanname'] = scanname
    parameters['scannumbers'] = scannumbers
    parameters['destpath'] = os.path.join(resultsdir,resultssubdir,out)
    if not instance.isarray(cfgfiles):
        cfgfiles = [cfgfiles]
    parameters['cfgfiles'] = [os.path.join(resultsdir,cfg) for cfg in cfgfiles]
    
    # Quantification
    parameters["qxrfgeometry"] = geometry
    
    # Image aligment
    if not fluoxas and not multi:
         parameters['alignmethod'] = None
    if instance.isstring(refimageindex):
        if refimageindex=='first':
            refimageindex = 0
        elif refimageindex=='middle':
            refimageindex = len(list(listtools.flatten(scannumbers)))/2
        elif refimageindex=='last':
            refimageindex = -1
        else:
            refimageindex = None # pair-wise
    elif instance.isnumber(refimageindex):
        if not instance.isinteger(refimageindex):
            refimageindex = max(0,min(refimageindex,1))
            refimageindex = int((len(list(listtools.flatten(scannumbers)))-1)*refimageindex)
    else:
        refimageindex = None # pair-wise
    parameters['refimageindex'] = refimageindex

    # Process
    h5filelast = process(**parameters)
    
    # Create EDF's for a single maps
    if not fluoxas:
        exportedf(h5filelast,**parameters)

def exportedf(h5name,**parameters):
    logger.info("EDF export {}:".format(h5name))
    
    instrument = create_hdf5_imagestacks.getinstrument(parameters)

    path = os.path.basename(h5name)
    n = 0
    while len(path)!=n:
        n = len(path)
        path = os.path.splitext(path)[0]
    path = os.path.join(os.path.dirname(h5name),"{}_results".format(path))
    
    # not necessary but clean in case of reruns
    if os.path.isdir(path):
        shutil.rmtree(path)
            
    filename = os.path.splitext(os.path.basename(h5name))[0]

    stacks, axes, procinfo = get_hdf5_imagestacks(h5name,instrument.h5stackgroups)
    counters = instrument.counters()

    ndet = sum(1 for k in stacks if k.name.startswith("detector"))

    with h5py.File(h5name) as hdf5FileObject:
        for g in stacks:
            if ndet==1:
                outpath = path
            else:
                outpath = os.path.join(path,g)
            
            if not os.path.exists(outpath):
                os.makedirs(outpath)
        
            for s in stacks[g]:
                if s in counters:
                    continue
                energy = hdf5FileObject[g][s][instrument.edfheaderkeys["energylabel"]]
                n = len(energy)
                for i in range(n):
                    outfile = s.split("/")[-1]
                    outfile = outfile.replace("-","_")
                    outfile = outfile.replace("(","")
                    outfile = outfile.replace(")","")
                    if n==1:
                        outfile = os.path.join(outpath,"{}.edf".format(outfile))
                    else:
                        outfile = os.path.join(outpath,"{}_{}{}.edf".format(outfile,energy[i],instrument.edfheaderkeys["energyunit"]))

                    logger.info(outfile)
                    edf.saveedf(outfile,np.squeeze(hdf5FileObject[g][s]["data"][...,i]),{'Title': s},overwrite=True)
                    

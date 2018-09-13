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

import json
import numpy as np
import logging
import re
import shutil
import os
import collections

from ..io import xiaedf
from ..io import edf
from ..io import spec
from ..io import nexus
from ..xrf.fit import PerformBatchFit
from ..common import units
from ..common import listtools

logger = logging.getLogger(__name__)



def detectorname(detector):
    # detector: "00", "S0"
    if detector:
        if detector.isdigit():
            name = "detector{:d}".format(int(detector))
        elif detector.startswith("S"):
            name = "detectorsum"
        else:
            raise "Unexpected detector name {}".format(detector)
    else:
        name = "counters"
    return name

def transmission_func(fluxt,flux0):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.divide(fluxt,flux0)

def absorbance_func(transmission):
    with np.errstate(divide='ignore', invalid='ignore'):
        return -np.log(np.clip(transmission,0,1))

def xrfnorm_func(xrf,flux,fluxref,xiaimage,detnr):
    if fluxref:
        norm = fluxref/xiaedf.normalizer(flux)
    else:
        norm = 1
    if xiaimage:
        xiaimage.onlyicrocr(True)
        stats = xiaimage.stats
        dtcor = xiaedf.deadtimecorrector(stats[...,0,detnr],stats[...,1,detnr])
        dtcor = dtcor.reshape(xrf.shape)
    else:
        dtcor = 1
    
    return xrf*norm*dtcor

def nanaverage(x):
    return np.nanmean(list(x),axis=0)

def nansum(x):
    return np.nansum(list(x),axis=0)

def nanmax(x):
    return np.nanmax(list(x),axis=0)

def identity(x):
    return x

class LazyArgument(object):

    def __init__(self,arg):
        self._arg = arg

    def data(self,*args):
        return self._arg 
    
    def __repr__(self):
        return self._arg
    
    def __str__(self):
        return self.__repr__()
        
class LazyArgumentEdf(LazyArgument):

    def __init__(self,filename):
        self._filename = filename

    def data(self,*args):
        return edf.edfimage(self._filename).data
    
    def __repr__(self):
        return self._filename
    
class LazyArgumentH5Dataset(LazyArgument):

    def __init__(self,groups):
        self._groups = groups

    def data(self,root,islice,stackdim):
        grp = root
        for g in self._groups:
            grp = grp[g]
        
        dset = grp[grp.attrs["signal"]]
        
        if stackdim == 0:
            data = dset[islice,...]
        elif stackdim == 1:
            data = dset[:,islice,:]
        else:
            data = dset[...,islice]
                            
        return data

    def __repr__(self):
        return '"{}"'.format('/'.join(self._groups))
    
    def __str__(self):
        return self.__repr__()

class LazyStackSlice(LazyArgument):

    def __init__(self,func=None,unpackargs=True):
        if func is None:
            self._func = identity
        else:
            self._func = func
        self._args = []
        self._unpackargs = unpackargs
        
    def data(self,*info):
        if self._unpackargs:
            return self._func(*list(self._arggen(*info)))
        else:
            return self._func(self._arggen(*info))
        
    def _arggen(self,*info):
        for x in self._args:
            if isinstance(x,LazyArgument):
                yield x.data(*info)
            else:
                yield x
                
    def appendarg(self,arg):
        self._args.append(arg)
           
    def appendarg_edf(self,filename):
        self.appendarg(LazyArgumentEdf(filename))

    def appendarg_h5dataset(self,*groups):
        self.appendarg(LazyArgumentH5Dataset(groups))

    def __repr__(self):
        if hasattr(self._func,"__name__"):
            func = self._func.__name__
        else:
            func = self._func
    
        return "{}({})".format(func,",".join([str(arg) for arg in self._args]))

class ImageStack(object):

    def __init__(self,jsonfile,qxrfgeometry=None):
        self.jsonfile = jsonfile
        self.qxrfgeometry = qxrfgeometry
    
    def create(self):
        # Processing configuration
        with open(self.jsonfile,'r') as f:
            self.config = json.load(f)

        self._preparestacks()
        self._processstacks()

    def _preparestacks(self):
        self._prepareinput()
        self._prepareoutput()
        self._prepare_processing_general()

        # Add stacks
        self._stack_add_counters()
        self._prepare_processing_xrfquant()
        self._stack_add_flux()
        self._correct_xrf_spectra()
        self._stack_add_xrffit()
        
        # Post processing
        self._postcorrect()
        self._sort_stacks()
    
    def _processstacks(self):
        with nexus.File(self.config["hdf5output"],mode='w') as f:
            # Save stack axes values
            self.axes = nexus.createaxes(f,self.stackaxes)

            # Save groups
            self._exportgroups(f)

            # Save stackinfo
            coordgrp = nexus.newNXentry(f,"stackinfo")
            for k in self.stackinfo:
                coordgrp[k] = self.stackinfo[k]

            # Add processing info
            nexus.addinfogroup(f,"fromraw",self.procinfo)

        logger.info("Saved {}".format(self.config["hdf5output"]))
        
    def _prepareinput(self):
        # Check data
        npaths = len(self.config["sourcepath"])
        if npaths != len(self.config["scanname"]):
            raise ValueError("Number of scan names must be the same as number of source paths.")
        if npaths != len(self.config["scannumbers"]):
            raise ValueError("Number of scan numbers must be the same as number of source paths.")
        
        # Get image stack
        self.xiastackraw = xiaedf.xiastack_mapnumbers(self.config["sourcepath"],self.config["scanname"],self.config["scannumbers"])
        if self.xiastackraw.isempty:
            raise IOError("Cannot find data: {}".format(self.xiastackraw.filedescription))

        # Exclude detectors
        ndetorg = self.xiastackraw.dshape[-1]
        self.xiastackraw.skipdetectors(self.config["exclude_detectors"])
        self.xiastackraw.keepdetectors(self.config["include_detectors"])
        includeddetectors = self.xiastackraw.xiadetectorselect_numbers(range(ndetorg))
        self.allowedgroups = ['detector{}'.format(i) for i in includeddetectors]+['counters','detectorsum']
        self.detectorindex = {idet:i for i,idet in enumerate(includeddetectors)}
    
        # Counter directory relative to the XIA files
        self.xiastackraw.counter_reldir(self.config["counter_reldir"])

    @property
    def nstack(self):
        #nstack, nrow, ncol, nchan, ndet = self.xiastackraw.dshape
        return self.xiastackraw.dshape[0]
    
    @property
    def ndet(self):
        return self.xiastackraw.dshape[-1]
        
    def _prepareoutput(self):
        self.stacks = collections.OrderedDict()
        self.stackaxes = [None]*3
        self.stackinfo = {}

        self.outstackdim = self.config["stackdim"]
        if self.outstackdim == 0:
            self.outimgdim = [1,2]
        elif self.outstackdim == 1:
            self.outimgdim = [0,2]
        else:
            self.outimgdim = [0,1]

        self.procinfo = {'config':self.jsonfile}
        
    def _stack_add_counters(self):
        # Check counters
        countersfound = set(self.xiastackraw.counterbasenames())
        counters = countersfound.intersection(self.config["counters"])

        if self.config["metadata"]=="xia":
            metacounters = "xia"
        else:
            if countersfound:
                metacounters = next(iter(countersfound))
            else:
                logger.warning("Metacounters for {} are not found".format(self.xiastackraw)) 
                metacounters = []
                
        # Extract metadata and counters from raw stack
        for imageindex,xiaimage in enumerate(self.xiastackraw):
            header = xiaimage.header(source=metacounters)
            motfast,motslow,sfast,sslow,time = self._getscanparameters(header)

            # Prepare axes and stackinfo
            if imageindex==0:
                for mot in self.config["stackinfo"]:
                    if mot != motfast and mot != motslow and mot in header:
                        self.stackinfo[mot] = np.full(self.nstack,np.nan)
                        
                self.stackaxes[self.outimgdim[1]] = sfast
                self.stackaxes[self.outimgdim[0]] = sslow
                self.stackaxes[self.outstackdim] = {"name":str(self.config["stacklabel"]),"data":np.full(self.nstack,np.nan,dtype=np.float32)}
                self.stackinfo["expotime"] = np.full(self.nstack,np.nan,dtype=np.float32)
                if self.fluxnorm:
                    self.stackinfo["xrfdetectorposition"] = np.full((self.nstack,self.ndetfit),np.nan,dtype=np.float32)

            # Add stackinfo
            for mot in self.stackinfo:
                if mot in header:
                    self.stackinfo[mot][imageindex] = np.float(header[mot])
                              
            # Add stack value
            if self.config["stacklabel"] in header:
                self.stackaxes[self.outstackdim]["data"][imageindex] = np.float(header[self.config["stacklabel"]])
            else:
                logger.warning("No stack counter in header (set to NaN)")
            
            # Add time
            self.stackinfo["expotime"][imageindex] = time
            
            # Lazy add counters
            files = xiaimage.ctrfilenames_used(counters)
            files = xiaedf.xiagroupdetectors(files)
            for detector,v1 in files.items():
                name = detectorname(detector)
                if name not in self.allowedgroups:
                    continue
                    
                # Prepare list of files
                if imageindex==0:
                    self.stacks[name] = collections.OrderedDict()
                    for ctr in v1:
                        self.stacks[name][ctr] = [None]*self.nstack
                        
                for ctr,f in v1.items():
                    # Add counter file
                    o = LazyStackSlice()
                    o.appendarg_edf(f[0])
                    self.stacks[name][ctr][imageindex] = o

    def _getscanparameters(self,header):
        """Get scan dimensions from header
        """
        result = {"name":"unknown"}
        
        if "speccmdlabel" in self.config:
            if self.config["speccmdlabel"] in header:
                o = spec.cmd_parser()
                result = o.parse(header[self.config["speccmdlabel"]])
        elif "fastlabel" in config and "slowlabel" in self.config:
            o = spec.edfheader_parser(fastlabel=self.config["fastlabel"],slowlabel=self.config["slowlabel"])
            result = o.parse(header)

        if result["name"]=="zapimage":
            sfast = {"name":result["motfast"],"data":spec.zapline_values(result["startfast"],result["endfast"],result["npixelsfast"])} 
            sslow = {"name":result["motslow"],"data":spec.ascan_values(result["startslow"],result["endslow"],result["nstepsslow"])} 
            if "time" in result:
                time = units.umagnitude(result["time"],units="s")
            else:
                time = np.nan
        else:
            logger.warning("No motor positions in header (using pixels).")
            sfast = {"name":"fast","data":np.arange(int(header["Dim_2"]))}
            sslow = {"name":"slow","data":np.arange(int(header["Dim_1"]))}
            time = np.nan
            
        return sfast["name"],sslow["name"],sfast,sslow,time
    
    def _prepare_processing_general(self):
        self.adddetectors = self.ndet>1 and self.config["adddetectors"]
        self.addspectra = self.adddetectors and self.config["addbeforefit"]
        if self.addspectra:
            self.ndetfit = 1
        else:
            self.ndetfit = self.ndet
            
        self.fluxnorm = self.qxrfgeometry is not None
        self.fluxnormbefore = self.fluxnorm and self.config["correctspectra"]
        self.fluxnormafter = self.fluxnorm and not self.fluxnormbefore
        self.procinfo["fluxnorm"] = self.fluxnorm
        
        dtcor = self.config["dtcor"] 
        self.dtcorbefore = dtcor and (self.config["correctspectra"] or self.addspectra)
        self.dtcorafter = dtcor and not self.dtcorbefore
        self.procinfo["dtneeded"] = False

        self.fit = self.config["fit"]

        self.group_detadd_functions = {}
        self.group_xrfcorrect = set()
        
    def _prepare_processing_xrfquant(self):
        if not self.fluxnorm:
            return 
            
        ngeometries = len(self.qxrfgeometry.xrfgeometries)
        if ngeometries!=self.ndetfit:
            raise RuntimeError("You need {} detector geometries, {} provides.".format(self.ndetfit,ngeometries))
        
        tile = lambda x: np.tile(np.array(x)[np.newaxis,:],(self.nstack,1))

        self.stackinfo["refflux"] = np.full(self.nstack,np.nan,dtype=np.float32)
        self.stackinfo["refexpotime"] = np.full(self.nstack,np.nan,dtype=np.float32)
        self.stackinfo["activearea"] = tile(self.qxrfgeometry.getxrfactivearea())
        self.stackinfo["anglein"] = tile(self.qxrfgeometry.getxrfanglein())
        self.stackinfo["angleout"] = tile(self.qxrfgeometry.getxrfangleout())
        self.stackinfo["sampledetdistance"] = np.full((self.nstack,self.ndetfit),np.nan,dtype=np.float32)
        if self.fluxnormbefore:
            self.stackinfo["I0toflux"] = [None]*self.nstack
        
        for imageindex,xiaimage in enumerate(self.xiastackraw):
            energy = self.stackaxes[self.outstackdim]["data"][imageindex]
            if not np.isnan(energy):
                time = self.stackinfo["expotime"][imageindex]
                if np.isnan(time):
                    time = None
                    
                xrfnormop,\
                self.stackinfo["refflux"][imageindex],\
                self.stackinfo["refexpotime"][imageindex],\
                self.stackinfo["expotime"][imageindex]\
                 =self.qxrfgeometry.xrfnormop(energy,expotime=time)
                
                if self.fluxnormbefore:
                    xiaimage.localnorm(self.config["fluxcounter"],func=xrfnormop)
                    self.stackinfo["I0toflux"][imageindex] = str(xrfnormop)
                    
                pos = self.stackinfo["xrfdetectorposition"][imageindex,:]
                if np.isfinite(pos).all():
                    self.qxrfgeometry.setxrfposition(pos)
                self.stackinfo["xrfdetectorposition"][imageindex,:] = self.qxrfgeometry.getxrfdetectorposition()
                self.stackinfo["sampledetdistance"][imageindex,:] = self.qxrfgeometry.getxrfdistance()
        
    def _correct_xrf_spectra(self):
        self.xiastackraw.detectorsum(self.addspectra)
        self.xiastackraw.dtcor(self.dtcorbefore)
        if self.dtcorbefore or self.addspectra or self.fluxnormbefore:
            label = ""
            if self.dtcorbefore:
                label = label+"dt"
            if self.fluxnormbefore:
                label = label+"fl"
            if label:
                label = label+"cor"
                radix = ["{}_{}".format(radix,label) for radix in self.config["scanname"]]
            else:
                radix = self.config["scanname"]
            
            # not necesarry but clean in case of re-runs
            if os.path.isdir(self.config["outdatapath"]):
                shutil.rmtree(self.config["outdatapath"])
            self.xiastackproc = xiaedf.xiastack_mapnumbers(self.config["outdatapath"],radix,self.config["scannumbers"])
            self.xiastackproc.overwrite(True)
            
            if self.addspectra:
                xialabels = ["xiaS1"]
            else:
                xialabels = self.xiastackraw.xialabels_used

            logger.info("Creating corrected XRF spectra ...")
            #xiastackproc.save(self.xiastackraw.data,xialabels=xialabels) # all in memory
            self.xiastackproc.save(self.xiastackraw,xialabels=xialabels) # least memory usage possible
        else:
            self.xiastackproc = self.xiastackraw
            
    def _stack_add_flux(self): 
        if not self.fluxnorm or ("fluxcounter" not in self.config and "transmissioncounter" not in self.config):
            return
            
        self.stackinfo["I0toflux"] = [None]*self.nstack
        self.stackinfo["Ittoflux"] = [None]*self.nstack        
        for imageindex in range(self.nstack):
            energy =self. stackaxes[self.outstackdim]["data"][imageindex] # handle missing energies specifically?
            name = detectorname(None)
            time = self.stackinfo["refexpotime"][imageindex]
            
            if "fluxcounter" in self.config:
                op,_ = self.qxrfgeometry.I0op(energy,expotime=time,removebeamfilters=False)
                if "calc_flux0" not in self.stacks[name]:
                    self.stacks[name]["calc_flux0"] = [None]*self.nstack
                o = LazyStackSlice(func=op)
                o.appendarg_h5dataset(name,self.config["fluxcounter"])
                self.stacks[name]["calc_flux0"][imageindex] = o
                self.stackinfo["I0toflux"][imageindex] = str(op)
                
            if "transmissioncounter" in self.config:
                op,_ = self.qxrfgeometry.Itop(energy,expotime=time,removebeamfilters=True)
                if "calc_fluxt" not in self.stacks[name]:
                    self.stacks[name]["calc_fluxt"] = [None]*self.nstack
                    self.stacks[name]["calc_transmission"] = [None]*self.nstack
                    self.stacks[name]["calc_absorbance"] = [None]*self.nstack
                
                o = LazyStackSlice(func=op)
                o.appendarg_h5dataset(name,self.config["transmissioncounter"])
                self.stacks[name]["calc_fluxt"][imageindex] = o
                self.stackinfo["Ittoflux"][imageindex] = str(op)     
                               
                o = LazyStackSlice(func=transmission_func)
                o.appendarg_h5dataset(name,"calc_fluxt")
                o.appendarg_h5dataset(name,"calc_flux0")
                self.stacks[name]["calc_transmission"][imageindex] = o
                
                o = LazyStackSlice(func=absorbance_func)
                o.appendarg_h5dataset(name,"calc_transmission")
                self.stacks[name]["calc_absorbance"][imageindex] = o
    
    def _stack_add_xrffit(self):
        # Fit data and add elemental maps
        if not self.fit:
            return
            
        logger.info("Fit XRF spectra ...")

        # not necessary but clean in case of re-runs
        if os.path.isdir(self.config["outfitpath"]):
            shutil.rmtree(self.config["outfitpath"])
                    
        if len(self.config["detectorcfg"])==1:
            fitcfg = self.config["detectorcfg"]*self.ndetfit
        else:
            fitcfg = self.config["detectorcfg"]
            if len(fitcfg)!=self.ndetfit:
                raise RuntimeError("You need {} configuration files, {} provides.".format(self.ndetfit,len(fitcfg)))
        
        for imageindex,xiaimage in enumerate(self.xiastackproc):
            if self.fluxnorm:
                quants = [{"time":self.stackinfo["refexpotime"][imageindex],\
                           "flux":self.stackinfo["refflux"][imageindex],\
                           "area":self.stackinfo["activearea"][imageindex,i],\
                           "anglein":self.stackinfo["anglein"][imageindex,i],\
                           "angleout":self.stackinfo["angleout"][imageindex,i],\
                           "distance":self.stackinfo["sampledetdistance"][imageindex,i]} for i in range(self.ndetfit)]
            else:
                quants = [{}]*self.ndetfit

            filestofit = xiaimage.datafilenames_used()
            filestofit = xiaedf.xiagroupdetectors(filestofit)
            for detector,cfg,quant in zip(filestofit,fitcfg,quants):
                # Fit
                outname = "{}_xia{}_{:04d}_0000".format(xiaimage.radix,detector,xiaimage.mapnum)
                energy = self.stackaxes[self.outstackdim]["data"][imageindex]

                files, labels = PerformBatchFit(filestofit[detector]["xia"],
                                       self.config["outfitpath"],outname,cfg,energy,
                                       fast=self.config["fastfitting"],mlines=self.config["mlines"],quant=quant)
                
                # Prepare list of files    
                name = detectorname(detector)
                if imageindex==0:
                    if name not in self.stacks:
                        self.stacks[name] = {}
                    for label in labels:
                        self.stacks[name][label] = [None]*self.nstack

                # Add file name           
                for f,label in zip(files,labels):
                    o = LazyStackSlice()
                    o.appendarg_edf(f)
                    self.stacks[name][label][imageindex] = o

    def _postcorrect(self):
        if self.dtcorafter:
            self._apply_postcorrect_xrf(self.fluxnormafter,self.dtcorafter)
            
        if self.adddetectors:
            self._add_detectors()

        if self.fluxnormafter and not self.dtcorafter:
            self._apply_postcorrect_xrf(self.fluxnormafter,self.dtcorafter)

    def _is_single_detector(self,grpname):
        return re.match("^detector([0-9]+)$",grpname)

    def _is_detector(self,grpname):
        return re.match("^detector([0-9]+|sum)$",grpname)

    def _add_detectors(self):
        detectors = [grpname for grpname in self.stacks if self._is_single_detector(grpname)]
        if "detectorsum" not in self.stacks:
            self.stacks["detectorsum"] = {}
        
        detectorsum_orggrps = list(self.stacks["detectorsum"].keys())
        
        for k1 in detectors:
            for k2 in self.stacks[k1]:
                if k2 in detectorsum_orggrps:
                    continue
                if k2 not in self.stacks["detectorsum"]:
                    func = self._group_add_detectors(k2)
                    self.stacks["detectorsum"][k2] = [LazyStackSlice(func=func,unpackargs=False) for _ in range(self.nstack)]
                for imageindex,arg in enumerate(self.stacks[k1][k2]):
                    self.stacks["detectorsum"][k2][imageindex].appendarg(arg)
            self.stacks.pop(k1)

    def _group_add_detectors(self,grpname):
        counter = any(k in grpname for k in self.config["counters"])
    
        if counter:
            func = sum
        elif grpname.startswith('w'):
            func = nanaverage
        elif 'chisq' in grpname:
            func = nanmax
        else:
            func = sum
        
        return func
        
    def _apply_postcorrect_xrf(self,fluxnorm,dtcor):
        detectors = [grpname for grpname in self.stacks if self._is_detector(grpname)]
        normname = detectorname(None)
        for k1 in detectors:
            for k2 in self.stacks[k1]:
                if not self._group_apply_postcorrect_xrf(k2):
                    continue

                keep = self.stacks[k1][k2]
                self.stacks[k1][k2] = [LazyStackSlice(func=xrfnorm_func) for _ in range(self.nstack)]

                for imageindex,(arg,xiaimage) in enumerate(zip(keep,self.xiastackproc)):
                    # arguments: xrf,flux,fluxref,xiaimage
                    self.stacks[k1][k2][imageindex].appendarg(arg)
                    if fluxnorm:
                        self.stacks[k1][k2][imageindex].appendarg_h5dataset(normname,"calc_flux0")
                        self.stacks[k1][k2][imageindex].appendarg(self.stackinfo["refflux"][imageindex])
                    else:
                        self.stacks[k1][k2][imageindex].appendarg(None)
                        self.stacks[k1][k2][imageindex].appendarg(None)
                    if dtcor and k1!="detectorsum":
                        self.stacks[k1][k2][imageindex].appendarg(xiaimage)
                        self.stacks[k1][k2][imageindex].appendarg(self.detectorindex[int(k1[8:])])
                    else:
                        self.stacks[k1][k2][imageindex].appendarg(None)
                        self.stacks[k1][k2][imageindex].appendarg(None)

    def _group_apply_postcorrect_xrf(self,grpname):
        return not any(k in grpname for k in self.config["counters"])
 
    def _sort_stacks(self):
        # Sort stack on stack axis value
        ind = np.argsort(self.stackaxes[self.outstackdim]["data"],kind='mergesort')
        self.stackaxes[self.outstackdim]["data"] = self.stackaxes[self.outstackdim]["data"][ind]
        for mot in self.stackinfo:
            try:
                self.stackinfo[mot] = self.stackinfo[mot][ind]
            except TypeError:
                self.stackinfo[mot] = listtools.listadvanced_int(self.stackinfo[mot],ind)

        for k1 in self.stacks:
            group = self.stacks[k1]
            for k2 in group:
                group[k2] = [group[k2][i] for i in ind]
    
    def _exportgroups(self,f):
        """Export groups of EDF stacks, summated or not
        """
        outshape = None
        for k1 in self.stacks: # detector or counter group
            if k1 in f:
                grp = f[k1]
            else:
                grp = nexus.newNXentry(f,k1)

            for k2 in self.stacks[k1]: # stack subgroup (Al-K, S-K, xmap_x1c, ...)
                nxdatagrp = None
                for imageindex,slicedata in enumerate(self.stacks[k1][k2]):
                    data = slicedata.data(f,imageindex,self.outstackdim)
                    
                    # Get destination for data
                    if k2 in grp:
                        dset = grp[k2][grp[k2].attrs["signal"]]
                    else:
                        logger.info("saving {}/{}".format(k1,k2))
                        nxdatagrp = nexus.newNXdata(grp,k2,"")
                        
                        # Allocate dataset
                        if outshape is None:
                            outshape = [0,0,0]
                            outshape[self.outimgdim[0]] = data.shape[0]
                            outshape[self.outimgdim[1]] = data.shape[1]
                            outshape[self.outstackdim] = len(self.stacks[k1][k2])
                        dset = nexus.createNXdataSignal(nxdatagrp,shape=outshape,chunks = True,dtype = np.float32)
                        dset[:] = np.nan

                        # Link axes to the new NXdata group
                        nexus.linkaxes(f,self.axes,[nxdatagrp])

                    logger.debug("Slice generator (index {}): {}".format(imageindex,slicedata))

                    # Some rows too much or rows missing:
                    if outshape[self.outimgdim[0]] > data.shape[0]:
                        data = np.pad(data,((0,outshape[self.outimgdim[0]]-data.shape[0]),(0,0)),'constant',constant_values=0)
                    elif outshape[self.outimgdim[0]] < data.shape[0]:
                        data = data[0:outshape[self.outimgdim[0]],:]

                    # Save data
                    if self.outstackdim == 0:
                        dset[imageindex,...] = data
                    elif self.outstackdim == 1:
                        dset[:,imageindex,:] = data
                    else:
                        dset[...,imageindex] = data

                else:
                    if nxdatagrp is not None:
                        # Replace subgroup k2 filename or calcinfo with NXentry name in stack
                        self.stacks[k1][k2] = nxdatagrp.name

        
def create_hdf5_imagestacks(jsonfile,qxrfgeometry=None):
    """Convert scanning data (XIA spectra + counters) to an HDF5 file:
        groups which contain NXdata classes
        3 axes datasets on the main level

    Returns:
        stacks(dict): {"counters": {"name1":lstack1,"name2":lstack2,...},
        "detector0":{"name3":lstack3,"name4":lstack4,...},
        "detector1":{"name3":lstack5,"name4":lstack6,...},...}
        lstack: an image stack given as an NXdata path

        axes(list): [{"name":"name1","fullname":"/axes/name1/data"},
        {"name":"name2","fullname":"/axes/name2/data"},
        {"name":"name3","fullname":"/axes/name3/data"}]
    """
    
    
    imagestack = ImageStack(jsonfile,qxrfgeometry=qxrfgeometry)
    imagestack.create()
    
    return imagestack.stacks,imagestack.axes,[imagestack.procinfo]


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

import numpy as np
import logging
import re
import shutil
import os
from collections import OrderedDict
import json

from ..io import xiaedf
from ..io import spec
from ..io import nexus
from ..io.utils import mkdir
from ..xrf.fit import PerformBatchFit
from ..common import units
from ..common import instance
from ..common import listtools
from ..instruments import configuration
from ..h5stacks.groups_hdf5_imagestacks import Group
from ..h5stacks import create_hdf5_lazy

logger = logging.getLogger(__name__)

def getinstrument(parameters):
    if "instrument" not in parameters:
        raise RuntimeError("You need to specify an instrument.")
    instrument = parameters["instrument"]
    if isinstance(instrument,configuration.InstrumentInfo):
        return instrument
    return configuration.factory(instrument,**parameters.get("instrument_parameters",{}))
    
def createconfig(sourcepath,destpath,scanname,scannumbers,cfgfiles,**parameters):
    
    instrument = getinstrument(parameters)
    stackdim = parameters.get("stackdim",2)
    dtcor = parameters.get("dtcor",True)
    fastfitting = parameters.get("fastfitting",True)
    adddetectors = parameters.get("adddetectors",True)
    addbeforefit = parameters.get("addbeforefit",True)
    mlines = parameters.get("mlines",{})
    exclude_detectors = parameters.get("exclude_detectors",None)
    include_detectors = parameters.get("include_detectors",None)
    noxia = parameters.get("noxia",False)
    encodercor = parameters.get("encodercor",{})
    qxrfgeometry = parameters.get("qxrfgeometry",None)
    qxrfgeometryparams = {"quantify":qxrfgeometry is not None}
    correctspectra = parameters.get("correctspectra",False)
    fluxid = parameters.get("fluxid","I0")
    transmissionid = parameters.get("transmissionid","It")
    #dtcorcounters = all(k in instrument.counterdict for k in ["xrficr","xrfocr"])
    
    if noxia:
        cfgfiles = None
    bfit = cfgfiles is not None

    if exclude_detectors is None:
        exclude_detectors = []

    if include_detectors is None:
        include_detectors = []

    if not instance.isarray(sourcepath):
        sourcepath = [sourcepath]
    if not instance.isarray(scanname):
        scanname = [scanname]
    if not instance.isarray(scannumbers):
        scannumbers = [scannumbers]
    if not instance.isarray(cfgfiles):
        cfgfiles = [cfgfiles]

    lst = []
    if noxia:
        lst.extend(["xrficr","xrfocr","xrfroi"])
    if not encodercor:
        lst.extend(["motors"])
    counters = instrument.counters(exclude=lst)
    counters.extend(parameters.get("counters",[]))
    
    config = {
            # Input
            "sourcepath": sourcepath,
            "counter_reldir": instrument.counter_reldir,
            "scanname": scanname,
            "scannumbers": scannumbers,
            "counters": counters,
            "fluxcounter": instrument.counterdict[fluxid+"_counts"],
            "transmissioncounter": instrument.counterdict[transmissionid+"_counts"],

            # Meta data
            "metadata": instrument.metadata,
            "stacklabel": instrument.edfheaderkeys["energylabel"],
            "speccmdlabel": instrument.edfheaderkeys["speclabel"],
            "fastlabel": instrument.edfheaderkeys["fastlabel"],
            "slowlabel": instrument.edfheaderkeys["slowlabel"],
            "timelabel": instrument.edfheaderkeys["timelabel"],
            "stackinfo": instrument.imagemotors,

            # Deadtime correction
            "dtcor": dtcor,
            "correctspectra": correctspectra,
            "adddetectors": adddetectors, # sum spectra
            "addbeforefit": addbeforefit, # sum fit results and detector counters
            "qxrfgeometry": qxrfgeometryparams,
            
            # Configuration for fitting
            "detectorcfg": cfgfiles,
            "mlines": mlines,
            "fit": bfit,
            "fastfitting": fastfitting,
            "exclude_detectors": exclude_detectors,
            "include_detectors": include_detectors,

            # Output directories
            "destpath": destpath,
            "outbase": scanname[0],
            "outdatapath": os.path.join(destpath,scanname[0]+"_data"),
            "outfitpath": os.path.join(destpath,scanname[0]+"_fit"),
            "hdf5output": os.path.join(destpath,scanname[0]+".h5"),
            "stackdim": stackdim
    }
    return config
    
class ImageStack(object):

    def __init__(self,sourcepath,destpath,scanname,scannumbers,cfgfiles,qxrfgeometry=None,**parameters):
        self.configin = {'sourcepath':sourcepath,
                         'destpath':destpath,
                         'scanname':scanname,
                         'scannumbers':scannumbers,
                         'cfgfiles':cfgfiles}
        self.configin.update(parameters)
        
        self.config = createconfig(sourcepath,destpath,scanname,scannumbers,cfgfiles,**parameters)
        self.qxrfgeometry = qxrfgeometry
    
    def create(self):
        self._preparestacks()
        self._processstacks()

    def _preparestacks(self):
        self._prepare_input()
        self._prepare_output()

        self._prepare_adddetector()
        self._prepare_fluxnorm()
        self._prepare_dtcor()
        self._stack_add_counters()
        self._prepare_xrfquant()
        self._prepare_xrffit()
        self._process_xiastackraw()
        self._stack_add_flux()
        self._stack_add_xrffit()

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
        
    def _prepare_input(self):
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

        dshape = self.xiastackraw.dshape #nstack, nrow, ncol, nchan, ndet 
        self.ndetorg = dshape[-1] # xia00, xia01, xiaS0, ...
        self.nstack = dshape[0]
        
    def _prepare_output(self):
        self.stacks = OrderedDict()
        self.stackaxes = [None]*3
        self.stackinfo = {}

        self.outstackdim = self.config["stackdim"]
        if self.outstackdim == 0:
            self.outimgdim = [1,2]
        elif self.outstackdim == 1:
            self.outimgdim = [0,2]
        else:
            self.outimgdim = [0,1]

        mkdir(self.config['destpath'])
        self._procinfo()
        
    def _procinfo(self):
        self.procinfo = {}
        destpath = self.config['destpath']
        outbase = self.config['outbase']
        
        jsonfile = os.path.join(destpath,outbase+".input.json")
        with open(jsonfile,'w') as f:
            json.dump(self.configin,f,indent=2)
        self.procinfo['input'] = jsonfile
        
        jsonfile = os.path.join(destpath,outbase+".process.json")
        with open(jsonfile,'w') as f:
            json.dump(self.config,f,indent=2)
        self.procinfo['process'] = jsonfile
        
        self.procinfo['result'] = {}
                
    def _prepare_adddetector(self):
        # Detector include/exclude
        self.xiastackraw.exclude_detectors(self.config["exclude_detectors"])
        include_detectors = list(listtools.flatten(self.config["include_detectors"]))
        self.xiastackraw.include_detectors(include_detectors)
        dshape = self.xiastackraw.dshape
        self.ndet = dshape[-1]

        # Do we need to add all detectors?
        adddetectors = (self.ndet>1) and self.config["adddetectors"]

        # Do we need to add detector groups?
        if any(len(instance.asarray(lst))>1 for lst in instance.asarray(self.config["include_detectors"])):
            self._include_detectors = {Group('S{:d}'.format(i)):instance.asarray(singledets)
                                       for i,singledets in enumerate(self.config["include_detectors"],1)}
        else:
            # These are not necessarily all detector!
            self._include_detectors = {Group('S1'):include_detectors}
        adddetectorgroups = len(self._include_detectors)>1

        # Do we need to add spectra (all or in groups)?
        self.addspectra = (adddetectors or adddetectorgroups) and self.config["addbeforefit"]
        self.xiastackraw.detectorsum(self.addspectra)
        
        # How many detectors need to be fitted?
        if self.addspectra:
            if adddetectorgroups:
                self.ndetfit = len(self._include_detectors)
            else:
                self.ndetfit = 1
        else:
            self.ndetfit = self.ndet

        # Sum after fitting:
        self._detectors_sumto = OrderedDict()
        if adddetectorgroups:
            # Add single detectors to their corresponding subgroup
            for dest,singledets in self._include_detectors.items():
                for det in singledets:
                    self._detectors_sumto[Group(det)] = dest
                
            if adddetectors:
                # Add subgroups to the sumgroup
                num = max(det.number for det in self._detectors_sumto)
                dest = Group('S{:d}'.format(num+1))
                for det in self._include_detectors:
                    self._detectors_sumto[Group(det)] = dest
        elif adddetectors:
            # Add single detectors to the sumgroup
            dest = Group('S1')
            for det in self.xiastackraw.xialabels_used:
                self._detectors_sumto[Group(det)] = dest

    def _prepare_dtcor(self):
        dtcor = self.config["dtcor"] 
        self.dtcorbefore = dtcor and (self.config["correctspectra"] or self.addspectra)
        self.dtcorafter = dtcor and not self.dtcorbefore
        self.procinfo["result"]["dtneeded"] = False
        self.xiastackraw.dtcor(self.dtcorbefore)

    def _prepare_fluxnorm(self):
        self.fluxnorm = self.qxrfgeometry is not None
        self.fluxnormbefore = self.fluxnorm and self.config["correctspectra"]
        self.fluxnormafter = self.fluxnorm and not self.fluxnormbefore
        self.procinfo["result"]["fluxnorm"] = self.fluxnorm

    def _prepare_xrfquant(self):  
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

    def _process_xiastackraw(self):
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
            
            # Only spectra, not counters (more efficient that way)
            if self.addspectra:
                logger.info("Creating corrected XRF spectra: {}".format(list(self._include_detectors.keys())))
                for det,include_detector in self._include_detectors.items():
                    self.xiastackraw.include_detectors(include_detector)
                    self.xiastackproc.save(self.xiastackraw,xialabels=[det.xialabel])
            else:
                xialabels = self.xiastackraw.xialabels_used
                logger.info("Creating corrected XRF spectra: {}".format(xialabels))
                self.xiastackproc.save(self.xiastackraw,xialabels=xialabels)
        else:
            logger.debug('Corrections (if any) are not applied to the XRF spectra')
            self.xiastackproc = self.xiastackraw
            
    def _stack_add_counters(self):
        self.xiastackraw.exclude_detectors(self.config["exclude_detectors"])
        self.xiastackraw.include_detectors(list(listtools.flatten(self.config["include_detectors"])))
        
        # Counter directory relative to the XIA files
        self.xiastackraw.counter_reldir(self.config["counter_reldir"])
        
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
        
        logger.info("Processing counters: {}".format(counters)) 
        
        # Extract metadata and counters from raw stack
        self.counters = set()
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
                name = Group(detector)

                # Prepare list of files
                if imageindex==0:
                    self.stacks[name] = OrderedDict()
                    for ctr in v1:
                        self.stacks[name][ctr] = [None]*self.nstack
                        
                for ctr,f in v1.items():
                    # Add counter file
                    o = create_hdf5_lazy.LazyStackSlice()
                    o.appendarg_edf(f[0])
                    self.stacks[name][ctr][imageindex] = o
                    self.counters.add(ctr)

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

    def _prepare_xrffit(self):
        self.fit = self.config["fit"]
        
    def _stack_add_flux(self): 
        if not self.fluxnorm or ("fluxcounter" not in self.config and "transmissioncounter" not in self.config):
            return
        
        self.stackinfo["I0toflux"] = [None]*self.nstack
        self.stackinfo["Ittoflux"] = [None]*self.nstack  
        name = Group(None)      
        for imageindex in range(self.nstack):
            energy = self.stackaxes[self.outstackdim]["data"][imageindex] # handle missing energies specifically?
            time = self.stackinfo["refexpotime"][imageindex]
            
            if "fluxcounter" in self.config:
                op,_ = self.qxrfgeometry.I0op(energy,expotime=time,removebeamfilters=False)
                if "calc_flux0" not in self.stacks[name]:
                    self.stacks[name]["calc_flux0"] = [None]*self.nstack
                o = create_hdf5_lazy.LazyStackSlice(func=op)
                o.appendarg_h5dataset(name,self.config["fluxcounter"])
                self.stacks[name]["calc_flux0"][imageindex] = o
                self.stackinfo["I0toflux"][imageindex] = str(op)
                
            if "transmissioncounter" in self.config:
                op,_ = self.qxrfgeometry.Itop(energy,expotime=time,removebeamfilters=True)
                if "calc_fluxt" not in self.stacks[name]:
                    self.stacks[name]["calc_fluxt"] = [None]*self.nstack
                    self.stacks[name]["calc_transmission"] = [None]*self.nstack
                    self.stacks[name]["calc_absorbance"] = [None]*self.nstack
                
                o = create_hdf5_lazy.LazyStackSlice(func=op)
                o.appendarg_h5dataset(name,self.config["transmissioncounter"])
                self.stacks[name]["calc_fluxt"][imageindex] = o
                self.stackinfo["Ittoflux"][imageindex] = str(op)     
                               
                o = create_hdf5_lazy.LazyStackSlice(func=create_hdf5_lazy.transmission_func)
                o.appendarg_h5dataset(name,"calc_fluxt")
                o.appendarg_h5dataset(name,"calc_flux0")
                self.stacks[name]["calc_transmission"][imageindex] = o
                
                o = create_hdf5_lazy.LazyStackSlice(func=create_hdf5_lazy.absorbance_func)
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
                logger.info('Pymca fit "detector{}" with "{}" ...'.format(detector,cfg))
                
                # Fit
                outname = "{}_xia{}_{:04d}_0000".format(xiaimage.radix,detector,xiaimage.mapnum)
                energy = self.stackaxes[self.outstackdim]["data"][imageindex]

                files, labels = PerformBatchFit(filestofit[detector]["xia"],
                                       self.config["outfitpath"],outname,cfg,energy,
                                       fast=self.config["fastfitting"],mlines=self.config["mlines"],quant=quant)
                
                # Prepare list of files    
                name = Group(detector)
                if imageindex==0:
                    if name not in self.stacks:
                        self.stacks[name] = {}
                    for label in labels:
                        self.stacks[name][label] = [None]*self.nstack

                # Add file name           
                for f,label in zip(files,labels):
                    o = create_hdf5_lazy.LazyStackSlice()
                    o.appendarg_edf(f)
                    self.stacks[name][label][imageindex] = o

    def _postcorrect(self):
        if self._detectors_sumto:
            self._apply_postcorrect_xrf(self.fluxnormafter,self.dtcorafter)
            self._add_detectors()
            
            #if self.dtcorafter:
            #    self._apply_postcorrect_xrf(False,self.dtcorafter)
            #self._add_detectors()
            #if self.fluxnormafter:
            #    self._apply_postcorrect_xrf(self.fluxnormafter,False)
        else:
            self._apply_postcorrect_xrf(self.fluxnormafter,self.dtcorafter)
            
    def _add_detectors(self):
        # Datasets in sum groups that are already their shouldn't be changed:
        fixedgroups = {}

        # Add detector to sum and remove it:
        for k1,sumdest in self._detectors_sumto.items():
            if k1 not in self.stacks:
                continue
                
            if sumdest not in self.stacks:
                self.stacks[sumdest] = {}
            if sumdest not in fixedgroups:
                fixedgroups[sumdest] = list(self.stacks[sumdest].keys())
                
            logger.debug('Remove {} and add to {}'.format(k1,sumdest))
            
            for k2 in self.stacks[k1]:
                if k2 in fixedgroups[sumdest]:
                    logger.debug('Do not add {}["{}"] to {}'.format(k1,k2,sumdest))
                    continue
                if k2 not in self.stacks[sumdest]:
                    func = self._group_add_detectors(k2)
                    self.stacks[sumdest][k2] = [create_hdf5_lazy.LazyStackSlice(func=func,unpackargs=False) for _ in range(self.nstack)]
                for imageindex,arg in enumerate(self.stacks[k1][k2]):
                    self.stacks[sumdest][k2][imageindex].appendarg(arg)
                    
            self.stacks.pop(k1)
        
    def _group_add_detectors(self,grpname):
        if grpname in self.counters:
            func = create_hdf5_lazy.sum_func
        elif grpname.startswith('w'):
            # mass fractions
            func = create_hdf5_lazy.nanmean_func
        elif 'chisq' in grpname:
            func = create_hdf5_lazy.nanmax_func
        else:
            # peak areas
            func = create_hdf5_lazy.sum_func
        
        return func
        
    def _apply_postcorrect_xrf(self,fluxnorm,dtcor):
        detectors = [detector for detector in self.stacks if detector.isdetector]
        normname = Group(None)

        logger.debug('Correction after XRF fitting (flux:{},dt:{})'.format(fluxnorm,dtcor))

        for k1 in detectors:
            for k2 in self.stacks[k1]:
                if k2 in self.counters:
                    continue
                
                keep = self.stacks[k1][k2]
                self.stacks[k1][k2] = [create_hdf5_lazy.LazyStackSlice(func=create_hdf5_lazy.xrfnorm_func) for _ in range(self.nstack)]
                for imageindex,(arg,xiaimage) in enumerate(zip(keep,self.xiastackproc)):
                    # arguments: xrf,flux,fluxref,xiaimage
                    self.stacks[k1][k2][imageindex].appendarg(arg)
                    if fluxnorm:
                        self.stacks[k1][k2][imageindex].appendarg_h5dataset(normname,"calc_flux0")
                        self.stacks[k1][k2][imageindex].appendarg(self.stackinfo["refflux"][imageindex])
                    else:
                        self.stacks[k1][k2][imageindex].appendarg(None)
                        self.stacks[k1][k2][imageindex].appendarg(None)
                    if dtcor and not k1.issum:
                        self.stacks[k1][k2][imageindex].appendarg(xiaimage)
                        self.stacks[k1][k2][imageindex].appendarg(k1.number)
                    else:
                        self.stacks[k1][k2][imageindex].appendarg(None)
                        self.stacks[k1][k2][imageindex].appendarg(None)
 
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
                    if k2 not in grp:
                        logger.info("saving {}/{}".format(k1,k2))
                    logger.debug("Slice generator (index {}): {}".format(imageindex,slicedata))
                    
                    data = slicedata.data(f,imageindex,self.outstackdim)
                    
                    # Get destination for data
                    if k2 in grp:
                        dset = grp[k2][grp[k2].attrs["signal"]]
                    else:
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
        
def create_hdf5_imagestacks(sourcepath,destpath,scanname,scannumbers,cfgfiles,**parameters):
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

    imagestack = ImageStack(sourcepath,destpath,scanname,scannumbers,cfgfiles,**parameters)
    imagestack.create()
    
    return imagestack.stacks,imagestack.axes,[imagestack.procinfo]


# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
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
from collections import OrderedDict
import logging

from . import nxtask
from . import nxresult
from . import nxlazy
from . import axis
from ..io import spec
from ..io import xiaedf
from ..io import localfs
from ..utils import listtools
from ..utils import instance
from ..utils import units
from ..xrf.fit import PerformBatchFit

logger = logging.getLogger(__name__)

class Task(nxtask.Task):
    
    DEFAULT_STACKDIM = 0
    
    def _parameters_defaults(self):
        super(Task,self)._parameters_defaults()
        self.allparams = [
                        # Input
                        "sourcepath",
                        "counter_reldir",
                        "scanname",
                        "scannumbers",
                        "counters",
                        "fluxcounter",
                        "transmissioncounter",
                        # Meta data
                        "metadata",
                        "stacklabel",
                        "speccmdlabel",
                        "fastlabel",
                        "slowlabel",
                        "timelabel",
                        "positionmotors",
                        "units",
                        # Data correction
                        "dtcor",
                        "correctspectra",
                        "adddetectors",
                        "addbeforefit",
                        "qxrfgeometry",
                        # Configuration for fitting
                        "detectorcfg",
                        "mlines",
                        "fit",
                        "fastfitting",
                        "exclude_detectors",
                        "include_detectors",
                        # Output directories (non hdf5 output)
                        "destpath",
                        "outbase",
                        "outdatapath",
                        "outfitpath",
                        ]
        self._required_parameters(*self.allparams)
        self.parameters['stackdim'] = self.parameters.get('stackdim',self.DEFAULT_STACKDIM)
        
        # TODO: temporary measure until pickleable
        self.qxrfgeometry = self.parameters.pop('qxrfgeometry')
        self.units = self.parameters.pop('units')

    def _parameters_filter(self):
        return super(Task,self)._parameters_filter()+self.allparams+['stackdim']

    def _execute(self):
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
        # NXdata positioners   
        positioners = self.nxresults.positioners()
        for ax in self.axes.values():
            positioners.add_axis(ax.name,ax.values,title=ax.title)

        # Processing info
        info = self.nxresults.nxcollection('info')
        for k,v in self.procinfo.items():
            info[k].mkfile(data=v)

        # Processing axes
        positioners = self.nxresults.nxcollection('stackaxes')
        for ax in self.infoaxes.values():
            positioners.add_axis(ax.name,ax.values,title=ax.title)

        # Processing datasets
        self._exportgroups()

    def _prepare_input(self):
        # Check data
        npaths = len(self.parameters["sourcepath"])
        if npaths != len(self.parameters["scanname"]):
            raise ValueError("Number of scan names must be the same as number of source paths.")
        if npaths != len(self.parameters["scannumbers"]):
            raise ValueError("Number of scan numbers must be the same as number of source paths.")
        
        # Get image stack
        self.xiastackraw = xiaedf.xiastack_mapnumbers(self.parameters["sourcepath"],self.parameters["scanname"],self.parameters["scannumbers"])
        if self.xiastackraw.isempty:
            raise IOError("Cannot find data: {}".format(self.xiastackraw.filedescription))

        dshape = self.xiastackraw.dshape #nstack, nrow, ncol, nchan, ndet 
        self.ndetorg = dshape[-1] # xia00, xia01, xiaS0, ...
        self.nstack = dshape[0]
        
    def _prepare_output(self):
        self.stacks = OrderedDict()
        self.axes_names = [None]*3
        self.axes = {}
        self.procinfo = {}
        self.infoaxes = {}
        
        self.outstackdim = self.parameters["stackdim"]
        if self.outstackdim == 0:
            self.outimgdim = [1,2]
        elif self.outstackdim == 1:
            self.outimgdim = [0,2]
        else:
            self.outimgdim = [0,1]
        
    def _prepare_adddetector(self):
        # Detector include/exclude
        self.xiastackraw.exclude_detectors = self.parameters["exclude_detectors"]
        include_detectors = list(listtools.flatten(self.parameters["include_detectors"]))
        self.xiastackraw.include_detectors = include_detectors
        dshape = self.xiastackraw.dshape
        self.ndet = dshape[-1]

        # Do we need to add all detectors?
        adddetectors = (self.ndet>1) and self.parameters["adddetectors"]

        # Do we need to add detector groups?
        if any(len(instance.asarray(lst))>1 for lst in instance.asarray(self.parameters["include_detectors"])):
            self._include_detectors = {nxresult.Group('S{:d}'.format(i)):instance.asarray(singledets)
                                       for i,singledets in enumerate(self.parameters["include_detectors"],1)}
        else:
            # These are not necessarily all detector!
            self._include_detectors = {nxresult.Group('S1'):include_detectors}
        adddetectorgroups = len(self._include_detectors)>1

        # Do we need to add spectra (all or in groups)?
        self.addspectra = (adddetectors or adddetectorgroups) and self.parameters["addbeforefit"]
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
                    self._detectors_sumto[nxresult.Group(det)] = dest
                
            if adddetectors:
                # Add subgroups to the sumgroup
                num = max(det.number for det in self._detectors_sumto)
                dest = nxresult.Group('S{:d}'.format(num+1))
                for det in self._include_detectors:
                    self._detectors_sumto[nxresult.Group(det)] = dest
        elif adddetectors:
            # Add single detectors to the sumgroup
            dest = nxresult.Group('S1')
            for det in self.xiastackraw.xialabels_used:
                self._detectors_sumto[nxresult.Group(det)] = dest

    def _prepare_dtcor(self):
        dtcor = self.parameters["dtcor"] 
        self.dtcorbefore = dtcor and (self.parameters["correctspectra"] or self.addspectra)
        self.dtcorafter = dtcor and not self.dtcorbefore
        self.procinfo["dtneeded"] = False
        self.xiastackraw.dtcor(self.dtcorbefore)

    def _prepare_fluxnorm(self):
        self.fluxnorm = self.qxrfgeometry is not None
        self.fluxnormbefore = self.fluxnorm and self.parameters["correctspectra"]
        self.fluxnormafter = self.fluxnorm and not self.fluxnormbefore
        self.procinfo["fluxnorm"] = self.fluxnorm

    def _prepare_xrfquant(self):  
        if not self.fluxnorm:
            return 
        ngeometries = len(self.qxrfgeometry.xrfgeometries)
        if ngeometries!=self.ndetfit:
            raise RuntimeError("You need {} detector geometries, {} provides.".format(self.ndetfit,ngeometries))
        
        tile = lambda x: np.tile(np.array(x)[np.newaxis,:],(self.nstack,1))

        self._add_info_axis("refflux",defaultunits='Hz')
        self._add_info_axis("refexpotime",defaultunits='s')
        self._add_info_axis("activearea",values=tile(self.qxrfgeometry.getxrfactivearea()),defaultunits='cm**2')
        self._add_info_axis("anglein",values=tile(self.qxrfgeometry.getxrfanglein()),defaultunits='deg')
        self._add_info_axis("angleout",values=tile(self.qxrfgeometry.getxrfangleout()),defaultunits='deg')
        self._add_info_axis("sampledetdistance",values=np.full((self.nstack,self.ndetfit),np.nan,dtype=np.float32),defaultunits='cm')
        
        if self.fluxnormbefore:
            self._add_info_axis("i0_to_norm_offset")
            self._add_info_axis("i0_to_norm_factor")
        
        for imageindex,xiaimage in enumerate(self.xiastackraw):
            energy = self.axes[self.axes_names[self.outstackdim]][imageindex].magnitude
            if not np.isnan(energy):
                time = self.infoaxes["expotime"][imageindex].magnitude
                if np.isnan(time):
                    time = None
                    
                xrfnormop,\
                self.infoaxes["refflux"][imageindex],\
                self.infoaxes["refexpotime"][imageindex],\
                self.infoaxes["expotime"][imageindex]\
                 =self.qxrfgeometry.xrfnormop(energy,expotime=time)
                
                if self.fluxnormbefore:
                    xiaimage.localnorm(self.parameters["fluxcounter"],func=xrfnormop)
                    self.infoaxes["i0_to_norm_offset"][imageindex] = xrfnormop.b
                    self.infoaxes["i0_to_norm_factor"][imageindex] = xrfnormop.m
                    
                pos = self.infoaxes["xrfdetectorposition"][imageindex,:]
                if np.isfinite(pos).all():
                    self.qxrfgeometry.setxrfposition(pos)
                self.infoaxes["xrfdetectorposition"][imageindex,:] = self.qxrfgeometry.getxrfdetectorposition()
                self.infoaxes["sampledetdistance"][imageindex,:] = self.qxrfgeometry.getxrfdistance()

    def _process_xiastackraw(self):
        if self.dtcorbefore or self.addspectra or self.fluxnormbefore:
            label = ""
            if self.dtcorbefore:
                label = label+"dt"
            if self.fluxnormbefore:
                label = label+"fl"
            if label:
                label = label+"cor"
                radix = ["{}_{}".format(radix,label) for radix in self.parameters["scanname"]]
            else:
                radix = self.parameters["scanname"]
            
            outpath = localfs.Path(self.parameters["outdatapath"])
            self.xiastackproc = xiaedf.xiastack_mapnumbers(outpath.path,radix,self.parameters["scannumbers"])
            shape = self.xiastackproc.dshape
            if shape:
                create = shape[:-1]!=self.xiastackraw.dshape[:-1]
            else:
                create = True
            if create:
                # Only spectra, not counters (more efficient that way)
                self.xiastackproc.overwrite(True)
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
                logger.info('Corrected XRF spectra already exist')
        else:
            logger.info('Corrections (if any) are not applied to the XRF spectra')
            self.xiastackproc = self.xiastackraw
        
    def _stack_add_counters(self):
        self.xiastackraw.exclude_detectors = self.parameters["exclude_detectors"]
        self.xiastackraw.include_detectors = list(listtools.flatten(self.parameters["include_detectors"]))
        
        # Counter directory relative to the XIA files
        self.xiastackraw.counter_reldir(self.parameters["counter_reldir"])
        
        # Check counters
        countersfound = set(self.xiastackraw.counterbasenames())
        counters = countersfound.intersection(self.parameters["counters"])

        if self.parameters["metadata"]=="xia":
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
            parsedheader = self._getscanparameters(header)

            # Prepare axes and stackinfo
            if imageindex==0:
                axes = parsedheader['axes']
                axesnames = [ax.name for ax in axes]
                for mot in self.parameters["positionmotors"]:
                    if mot not in axesnames:
                        self._add_stack_axis(mot,self.units.get(mot,'dimensionless'))
                self._add_grid_axis(axes[0],index=self.outimgdim[0])
                self._add_grid_axis(axes[1],index=self.outimgdim[1])
                self._add_stack_axis(self.parameters['stacklabel'],
                                     parsedheader['stackvalue'].units,
                                     index=self.outstackdim)
                self._add_info_axis("expotime",
                                    defaultunits=parsedheader['time'].units)
                if self.fluxnorm:
                    values = np.full((self.nstack,self.ndetfit),np.nan,dtype=np.float32)
                    self._add_info_axis("xrfdetectorposition",values=values,defaultunits='cm')
            
            # Add stackinfo
            for mot in self.axes:
                if mot in header:
                    self.axes[mot][imageindex] = np.float(header[mot])
                              
            # Add stack value
            self.axes[self.axes_names[self.outstackdim]][imageindex] = parsedheader['stackvalue'].magnitude

            # Add time
            self.infoaxes["expotime"][imageindex] = parsedheader['time'].magnitude
            
            # Lazy add counters
            files = xiaimage.ctrfilenames_used(counters)
            files = xiaedf.xiagroupdetectors(files)
            for detector,v1 in files.items():
                name = nxresult.Group(detector)

                # Prepare list of files
                if imageindex==0:
                    self.stacks[name] = OrderedDict()
                    for ctr in v1:
                        self.stacks[name][ctr] = [None]*self.nstack
                        
                for ctr,f in v1.items():
                    # Add counter file
                    o = nxlazy.LazyStackSlice()
                    o.appendarg_edf(f[0])
                    self.stacks[name][ctr][imageindex] = o
                    self.counters.add(ctr)
    
    def _add_grid_axis(self,axis,index=None):
        self.axes[axis.name] = axis
        if index is not None:
            self.axes_names[index] = axis.name
    
    def _add_stack_axis(self,name,unit,index=None):
        values = np.full(self.nstack,np.nan,dtype=np.float32)
        values = units.Quantity(values,units=unit)
        ax = axis.Axis(values,name=name)
        self._add_grid_axis(ax,index=index)

    def _add_info_axis(self,name,values=None,defaultunits=None,type=None):
        if values is None:
            values = np.full(self.nstack,np.nan,dtype=np.float32)
        values = units.Quantity(values,units=self.units.get(name,defaultunits))
        self.infoaxes[name] = axis.Axis(values,type=type,name=name)
        
    def _getscanparameters(self,header):
        """Get scan dimensions from header
        """
        kwargs = {}
        kwargs['speclabel'] = self.parameters.get('speccmdlabel',None)
        kwargs['slowlabel'] = self.parameters.get('slowlabel',None)
        kwargs['fastlabel'] = self.parameters.get('fastlabel',None)
        kwargs['stackvalue'] = self.parameters.get('stacklabel',None)
        kwargs['time'] = self.parameters.get('timelabel',None)
        o = spec.edfheader_parser(units=self.units,**kwargs)
        return o.parse(header)

    def _prepare_xrffit(self):
        self.fit = self.parameters["fit"]
        
    def _stack_add_flux(self): 
        if not self.fluxnorm or ("fluxcounter" not in self.parameters and "transmissioncounter" not in self.parameters):
            return
        
        self._add_info_axis("i0_to_flux_offset",defaultunits='Hz')
        self._add_info_axis("i0_to_flux_factor",defaultunits='Hz')
        self._add_info_axis("it_to_flux_offset",defaultunits='Hz')
        self._add_info_axis("it_to_flux_factor",defaultunits='Hz')
        
        name = nxresult.Group(None)
        for imageindex in range(self.nstack):
            energy = self.axes[self.axes_names[self.outstackdim]][imageindex].magnitude
            time = self.infoaxes["refexpotime"][imageindex]
            
            if "fluxcounter" in self.parameters:
                op,_ = self.qxrfgeometry.I0op(energy,expotime=time,removebeamfilters=False)
                if "calc_flux0" not in self.stacks[name]:
                    self.stacks[name]["calc_flux0"] = [None]*self.nstack
                o = nxlazy.LazyStackSlice(func=op)
                o.appendarg_h5dataset(self.nxresults[str(name)][self.parameters["fluxcounter"]])
                self.stacks[name]["calc_flux0"][imageindex] = o
                self.infoaxes["i0_to_flux_offset"][imageindex] = op.b
                self.infoaxes["i0_to_flux_factor"][imageindex] = op.m
                
            if "transmissioncounter" in self.parameters:
                op,_ = self.qxrfgeometry.Itop(energy,expotime=time,removebeamfilters=True)
                if "calc_fluxt" not in self.stacks[name]:
                    self.stacks[name]["calc_fluxt"] = [None]*self.nstack
                    self.stacks[name]["calc_transmission"] = [None]*self.nstack
                    self.stacks[name]["calc_absorbance"] = [None]*self.nstack
                
                o = nxlazy.LazyStackSlice(func=op)
                o.appendarg_h5dataset(self.nxresults[str(name)][self.parameters["transmissioncounter"]])
                self.stacks[name]["calc_fluxt"][imageindex] = o
                self.infoaxes["it_to_flux_offset"][imageindex] = op.b
                self.infoaxes["it_to_flux_factor"][imageindex] = op.m
                               
                o = nxlazy.LazyStackSlice(func=nxlazy.transmission_func)
                o.appendarg_h5dataset(self.nxresults[str(name)]["calc_fluxt"])
                o.appendarg_h5dataset(self.nxresults[str(name)]["calc_flux0"])
                self.stacks[name]["calc_transmission"][imageindex] = o
                
                o = nxlazy.LazyStackSlice(func=nxlazy.absorbance_func)
                o.appendarg_h5dataset(self.nxresults[str(name)]["calc_transmission"])
                self.stacks[name]["calc_absorbance"][imageindex] = o
        
    def _stack_add_xrffit(self):
        # Fit data and add elemental maps
        if not self.fit:
            return
            
        logger.info("Fit XRF spectra ...")

        # not necessary but clean in case of re-runs
        outpath = localfs.Path(self.parameters["outfitpath"])
        outpath.remove(recursive=True)

        if len(self.parameters["detectorcfg"])==1:
            fitcfg = self.parameters["detectorcfg"]*self.ndetfit
        else:
            fitcfg = self.parameters["detectorcfg"]
            if len(fitcfg)!=self.ndetfit:
                raise RuntimeError("You need {} configuration files, {} provides.".format(self.ndetfit,len(fitcfg)))
        
        for imageindex,xiaimage in enumerate(self.xiastackproc):
            if self.fluxnorm:
                quants = [{"time":self.infoaxes["refexpotime"][imageindex].to('s').magnitude,\
                           "flux":self.infoaxes["refflux"][imageindex].to('Hz').magnitude,\
                           "area":self.infoaxes["activearea"][imageindex,i].to('cm**2').magnitude,\
                           "anglein":self.infoaxes["anglein"][imageindex,i].to('deg').magnitude,\
                           "angleout":self.infoaxes["angleout"][imageindex,i].to('deg').magnitude,\
                           "distance":self.infoaxes["sampledetdistance"][imageindex,i].to('cm').magnitude} for i in range(self.ndetfit)]
            else:
                quants = [{}]*self.ndetfit

            filestofit = xiaimage.datafilenames_used()
            filestofit = xiaedf.xiagroupdetectors(filestofit)
            for detector,cfg,quant in zip(filestofit,fitcfg,quants):
                logger.info('Pymca fit "detector{}" with "{}" ...'.format(detector,cfg))
                
                # Fit
                outname = "{}_xia{}_{:04d}_0000".format(xiaimage.radix,detector,xiaimage.mapnum)
                energy = self.axes[self.axes_names[self.outstackdim]][imageindex].magnitude

                files, labels = PerformBatchFit(filestofit[detector]["xia"],
                                       outpath.path,outname,cfg,energy,
                                       fast=self.parameters["fastfitting"],mlines=self.parameters["mlines"],quant=quant)
                
                # Prepare list of files    
                name = nxresult.Group(detector)
                if imageindex==0:
                    if name not in self.stacks:
                        self.stacks[name] = {}
                    for label in sorted(labels):
                        self.stacks[name][label] = [None]*self.nstack

                # Add file name           
                for f,label in zip(files,labels):
                    o = nxlazy.LazyStackSlice()
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
                    self.stacks[sumdest][k2] = [nxlazy.LazyStackSlice(func=func,unpackargs=False) for _ in range(self.nstack)]
                for imageindex,arg in enumerate(self.stacks[k1][k2]):
                    self.stacks[sumdest][k2][imageindex].appendarg(arg)
                    
            self.stacks.pop(k1)
        
    def _group_add_detectors(self,grpname):
        if grpname in self.counters:
            func = nxlazy.sum_func
        elif grpname.startswith('w'):
            # mass fractions
            func = nxlazy.nanmean_func
        elif 'chisq' in grpname:
            func = nxlazy.nanmax_func
        else:
            # peak areas
            func = nxlazy.sum_func
        
        return func
        
    def _apply_postcorrect_xrf(self,fluxnorm,dtcor):
        detectors = [detector for detector in self.stacks if detector.isdetector]
        normname = nxresult.Group(None)

        logger.debug('Correction after XRF fitting (flux:{},dt:{})'.format(fluxnorm,dtcor))

        for k1 in detectors:
            for k2 in self.stacks[k1]:
                if k2 in self.counters:
                    continue
                
                keep = self.stacks[k1][k2]
                self.stacks[k1][k2] = [nxlazy.LazyStackSlice(func=nxlazy.xrfnorm_func) for _ in range(self.nstack)]
                for imageindex,(arg,xiaimage) in enumerate(zip(keep,self.xiastackproc)):
                    # arguments: xrf,flux,fluxref,xiaimage
                    self.stacks[k1][k2][imageindex].appendarg(arg)
                    if fluxnorm:
                        self.stacks[k1][k2][imageindex].appendarg_h5dataset(self.nxresults[str(normname)]["calc_flux0"])
                        self.stacks[k1][k2][imageindex].appendarg(self.infoaxes["refflux"][imageindex])
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
        mot = self.axes_names[self.outstackdim]
        ind = np.argsort(self.axes[mot].magnitude,kind='mergesort')
        self.axes[mot].values = self.axes[mot][ind]

        for name in self.infoaxes:
            self.infoaxes[name].values = self.infoaxes[name][ind,...]

        for k1 in self.stacks:
            group = self.stacks[k1]
            for k2 in group:
                group[k2] = [group[k2][i] for i in ind]
    
    def _exportgroups(self):
        """Export groups of EDF stacks, summed or not
        """
        outshape = None
        for k1 in self.stacks: # detector or counter group
            nxdata = self.nxresults.nxdata(str(k1))
            for k2 in self.stacks[k1]: # stack subgroup (Al-K, S-K, xmap_x1c, ...)
                signal = nxdata[str(k2)]
                for imageindex,slicedata in enumerate(self.stacks[k1][k2]):
                    if not signal.exists:
                        logger.info("saving {}/{}".format(k1,k2))
                    logger.debug("Slice generator (index {}): {}".format(imageindex,slicedata))
                    
                    data = slicedata.data(imageindex,self.outstackdim)
                    
                    # Allocate destination
                    if not signal.exists:
                        if outshape is None:
                            outshape = [0,0,0]
                            outshape[self.outimgdim[0]] = data.shape[0]
                            outshape[self.outimgdim[1]] = data.shape[1]
                            outshape[self.outstackdim] = len(self.stacks[k1][k2])
                        nxdata.add_signal(name=signal.name,shape=outshape,chunks=True,dtype=np.float32)

                    # Some rows too much or rows missing:
                    if outshape[self.outimgdim[0]] > data.shape[0]:
                        data = np.pad(data,((0,outshape[self.outimgdim[0]]-data.shape[0]),(0,0)),'constant',constant_values=0)
                    elif outshape[self.outimgdim[0]] < data.shape[0]:
                        data = data[0:outshape[self.outimgdim[0]],:]

                    # Save data
                    with signal.open() as dset:
                        if self.outstackdim == 0:
                            dset[imageindex,...] = data
                        elif self.outstackdim == 1:
                            dset[:,imageindex,:] = data
                        else:
                            dset[...,imageindex] = data

                if signal.exists:
                    self.stacks[k1][k2] = signal
                    
            if list(nxdata.signals):
                nxdata.set_axes(*self.axes_names)
                

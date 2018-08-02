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

import collections
import os

from .. import ureg
from ..common.classfactory import with_metaclass
from ..common import instance
from ..common import listtools

class InstrumentInfo(with_metaclass(object)):

    def __init__(self,**info):
        self.imagemotors = info.get("imagemotors",[])
        self.imageaxes = info.get("imageaxes",("y","x"))
        
        self.units = collections.defaultdict(lambda x:ureg.dimensionless)
        self.units.update(info.get("units",{}))
        
        self.encoderresolution = info.get("encoderresolution",{})# steps/motor unit
        self.compensationmotors = info.get("compensationmotors",{})
        
        self.edfheaderkeys = {"speclabel":"title",\
                              "energylabel":"energy",\
                              "timelabel":"timelabel",\
                              "energyunit":"keV",\
                              "timeunit":"s",\
                              "fastlabel":"fast",\
                              "slowlabel":"slow"}
        self.edfheaderkeys.update(info.get("edfheaderkeys",{}))

        self.counterdict = {"I0_counts":"I0","It_counts":"It"}
        self.counterdict.update(info.get("counterdict",[]))
        self.counter_reldir = info.get("counter_reldir",".")
        self.metadata = info.get("metadata","counters")
        
        self.xraysource = info.get("xraysource","synchrotron")
        self._diodeI0 = info.get("diodeI0",{})
        self._diodeIt = info.get("diodeIt",{})
        self._optics = info.get("optics",{})
        
        self._specmotornames = info.get("specmotornames",{})
        self.speccounternames = info.get("speccounternames",{})
        
        self.h5stackgroups = info.get("h5stackgroups",["counters","^detector(\d+|sum)$"])

        self.datalocation = info.get("datalocation",{})
        if "xrf_dynamic" not in self.datalocation:
            self.datalocation["xrf_dynamic"] = {"subdir":""}
        if "xrf_static" not in self.datalocation:
            self.datalocation["xrf_static"] = {"subdir":""}

    def counters(self,include=None,exclude=None):
        if include is None:
            lst = self.counterdict.keys()
        else:
            lst = include
        if exclude is None:
            exclude = []
        lst = [k for k in lst if k not in exclude]
            
        ret = [self.counterdict[k] for k in lst if k in self.counterdict]
        return list(listtools.flatten(ret))

    @property
    def metadata(self):
        return self._metadata

    @metadata.setter
    def metadata(self,value):
        if value=="xia":
            self._metadata = "xia"
        else:
            self._metadata = "counters"

    @staticmethod
    def _devicename(device,motordict):
        name = ""
        for k,f in device.items():
            name = f(motordict[k])
            if name:
                return name
        if "default" in device:
            name = device["default"]
        return name
    
    def diodeI0(self,motordict):
        return self._devicename(self._diodeI0,motordict)
    
    def diodeIt(self,motordict):
        return self._devicename(self._diodeIt,motordict)
    
    def optics(self,motordict):
        return self._devicename(self._optics,motordict)
            
    def specmotornames(self,names):
        return [self._specmotornames[k] for k in names]
    
    def motornames(self):
        return self._specmotornames.keys()

    def xrfsubdir(self,type="dynamic"):
        return self.datalocation["xrf_{}".format(type)]["subdir"]

    def xrflocation(self,sample,dataset,type="dynamic"):
        # root/sample/sample_dataset/xrfsubdir/sample_dataset*.*
        xrfsubdir = self.xrfsubdir(type=type)
        radix = "_".join([sample,dataset])
        subdir = os.path.join(sample,radix,xrfsubdir)
        return radix,subdir


class ESRF_ID21_SXM(InstrumentInfo):
    aliases = ["sxm","id21"]
    
    def __init__(self,**info):
        info["imagemotors"] = info.get("imagemotors",["samy","sampy","samz","sampz"])
        info["imageaxes"] = info.get("imageaxes",("z","y"))
        info["compensationmotors"] = info.get("compensationmotors",\
                                    {"samy":["sampy"],\
                                    "samz":["sampz"],\
                                    "sampy":["samy"],\
                                    "sampz":["samz"]})
        info["encoderresolution"] = info.get("encoderresolution",\
                                    {"samy":52500,"samz":50000})
                                    
        info["edfheaderkeys"] = info.get("edfheaderkeys",\
                                    {"speclabel":"Title",\
                                    "energylabel":"DCM_Energy"})
                                    
        info["counterdict"] = info.get("counterdict",\
                                    {"I0_counts":"arr_iodet",\
                                    "It_counts":"arr_idet",\
                                    "If_counts":"arr_fdet",\
                                    "xrficr":"xmap_icr",\
                                    "xrfocr":"xmap_ocr",\
                                    "encoders":["arr_samy","arr_samz"],\
                                    "xrfroi":["xmap_x1","xmap_x1c","xmap_x2","xmap_x2c","xmap_x3","xmap_x3c"],\
                                    "calc":["arr_absorp1","arr_absorp2","arr_absorp3"],\
                                    "counters":["arr_","xmap_"]})
        
        info["diodeI0"] = info.get("diodeI0",\
                                  {"istopz":lambda x:"iodet2" if abs(x+20)<abs(x+1.3) else "",\
                                   "ioz":lambda x:"iodet1" if abs(x-7)<abs(x-23) else ""}
                                  )

        info["diodeIt"] = info.get("diodeIt",\
                                  {"detz":lambda x:"idet" if abs(x-10)<abs(x-21) else ""}
                                  )
                                 
        info["optics"] = info.get("optics",\
                                  {"zpz": lambda x: "KB" if abs(x-7)<abs(x-6.5) else ""}
                                 )
        
        info["specmotornames"] = info.get("specmotornames",\
                                          {"ioz":"diodeIoZ",\
                                           "istopz":"istopz",\
                                           "energy":"Energy MONO",\
                                           "zpz":"zpz",\
                                           "detz":"detz",\
                                           }
                                         )
        
        info["speccounternames"] = info.get("speccounternames",\
                                          {"I0_counts":"iodet",\
                                           "It_counts":"idet",\
                                           "energy":"energym",\
                                           "time":"seconds",\
                                           "It_photons":"photons"}
                                         )
        
        info["units"] = info.get("units",{})                           
        info["units"]["time"] = info["units"].get("time",ureg.seconds)
        info["units"]["energy"] = info["units"].get("energy",ureg.Unit("keV"))
        info["units"]["I0_counts"] = info["units"].get("I0_counts",ureg.dimensionless)
        info["units"]["It_counts"] = info["units"].get("It_counts",ureg.dimensionless)
        info["units"]["It_photons"] = info["units"].get("It_photons",ureg.dimensionless)
        info["units"]["samy"] = info["units"].get("samy",ureg.millimeter)
        info["units"]["samz"] = info["units"].get("samz",ureg.millimeter)
        info["units"]["samx"] = info["units"].get("samx",ureg.millimeter)
        info["units"]["sampz"] = info["units"].get("sampz",ureg.micrometer)
        info["units"]["sampy"] = info["units"].get("sampy",ureg.micrometer)
                       
        info["datalocation"] = info.get("datalocation",{"xrf_dynamic":"zap","xrf_static":"xrf"})
             
        super(ESRF_ID21_SXM,self).__init__(**info)
                

class ESRF_ID21_MICRODIFF(InstrumentInfo):
    aliases = ["microdiff"]
    
    def __init__(self,**info):
        info["imagemotors"] = info.get("imagemotors",["samh","samph","samv","sampv"])
        info["imageaxes"] = info.get("imageaxes",("v","h"))
        info["units"] = info.get("units",\
                                    {"samh":ureg.millimeter,\
                                    "samv":ureg.millimeter,\
                                    "samd":ureg.millimeter,\
                                    "samph":ureg.micrometer,\
                                    "sampv":ureg.micrometer})
        info["compensationmotors"] = info.get("compensationmotors",\
                                    {"samh":["samph"],\
                                    "samv":["sampv"],\
                                    "samph":["samh"],\
                                    "sampv":["samv"]})
                                    
        info["counterdict"] = info.get("counterdict",\
                                    {"I0_counts":"zap_iodet",\
                                    "It_counts":"zap_idet",\
                                    "xrficr":"xmap_icr",\
                                    "xrfocr":"xmap_ocr",\
                                    "xrfroi":["xmap_x1","xmap_x1c","xmap_x2","xmap_x2c","xmap_x3","xmap_x3c"]})
        
        info["datalocation"] = info.get("datalocation",{"xrf_dynamic":"zap","xrf_static":"xrf","xrd_dynamic":"zap","xrd_static":"xrd"})

        super(ESRF_ID21_MICRODIFF,self).__init__(**info)


class ESRF_ID16B(InstrumentInfo):
    aliases = ["id16b"]
    
    def __init__(self,**info):
        info["imagemotors"] = info.get("imagemotors",["sy","sz","sampy","sampz"])
        info["imageaxes"] = info.get("imageaxes",("z","y"))
        info["units"] = info.get("units",\
                                    {"sx":ureg.millimeter,\
                                    "sy":ureg.millimeter,\
                                    "sz":ureg.millimeter,\
                                    "sampy":ureg.micrometer,\
                                    "sampz":ureg.micrometer})
        info["compensationmotors"] = info.get("compensationmotors",\
                                    {"sy":["sampy"],\
                                    "sz":["sampz"],\
                                    "sampy":["sy"],\
                                    "sampz":["sz"]})
                                    
        info["diodeI0"] = info.get("diodeI0",\
                                  {"default":"id16b_IC",\
                                   "unknown":"id16b_I0"},\
                                  )
                                    
        info["edfheaderkeys"] = info.get("edfheaderkeys",\
                                    {"speclabel":"title",\
                                    "energylabel":"energy"})

        info["diodeIt"] = info.get("diodeIt",\
                                  {"default":"id16b_It"}\
                                  )

        info["optics"] = info.get("optics",\
                                  {"zpz": "KB"}
                                 )
        
        info["counterdict"] = info.get("counterdict",\
                                    {"I0_counts":"zap_p201_I0",\
                                    "It_counts":"zap_p201_It",\
                                    "IC_counts":"zap_p201_IC",\
                                    "counters":["zap_","xmap_"]})

        info["speccounternames"] = info.get("speccounternames",\
                                          {"IC_counts":"p201_IC",\
                                           "I0_counts":"p201_I0",\
                                           "It_counts":"p201_It",\
                                           "IC_current":"IC",\
                                           "I0_current":"I0",\
                                           "It_current":"It",\
                                           "energy":"energy",\
                                           "time":"Seconds",\
                                           "It_flux":"flux_It"}
                                         )

        info["units"] = info.get("units",{})                           
        info["units"]["time"] = info["units"].get("time",ureg.seconds)
        info["units"]["I0_counts"] = info["units"].get("I0_counts",ureg.dimensionless)
        info["units"]["It_counts"] = info["units"].get("It_counts",ureg.dimensionless)
        info["units"]["IC_counts"] = info["units"].get("IC_counts",ureg.dimensionless)
        info["units"]["energy"] = info["units"].get("energy",ureg.Unit("keV"))
        info["units"]["It_flux"] = info["units"].get("It_flux",ureg.Unit("Hz"))
        info["units"]["I0_current"] = info["units"].get("I0_current",ureg.Unit("pA"))
        info["units"]["It_current"] = info["units"].get("It_current",ureg.Unit("pA"))
        info["units"]["IC_current"] = info["units"].get("IC_current",ureg.Unit("pA"))

        info["metadata"] = info.get("metadata","xia")

        super(ESRF_ID16B,self).__init__(**info)

    def xrflocation(self,sample,dataset,type="dynamic"):
        # root/sample/xrfsubdir/dataset*.*
        xrfsubdir = self.xrfsubdir(type=type)
        if not dataset:
            dataset = sample
        radix = dataset
        subdir = os.path.join(sample,xrfsubdir)
        return radix,subdir


        
factory = InstrumentInfo.factory


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

        self.counterdict = {"I0":"I0","It":"It"}
        self.counterdict.update(info.get("counterdict",[]))
        self.counter_reldir = info.get("counter_reldir",".")
        self.metadata = info.get("metadata","counters")

        self.h5stackgroups = info.get("h5stackgroups",["counters","^detector(\d+|sum)$"])

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
            
            
class ESRF_ID21_SXM(InstrumentInfo):
    aliases = ["sxm","id21"]
    
    def __init__(self,**info):
        info["imagemotors"] = info.get("imagemotors",["samy","sampy","samz","sampz"])
        info["imageaxes"] = info.get("imageaxes",("z","y"))
        info["units"] = info.get("units",\
                                    {"samy":ureg.millimeter,\
                                    "samz":ureg.millimeter,\
                                    "samx":ureg.millimeter,\
                                    "sampy":ureg.micrometer,\
                                    "sampz":ureg.micrometer})
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
                                    {"I0":"arr_iodet",\
                                    "It":"arr_idet",\
                                    "If":"arr_fdet",\
                                    "xrficr":"xmap_icr",\
                                    "xrfocr":"xmap_ocr",\
                                    "encoders":["arr_samy","arr_samz"],\
                                    "xrfroi":["xmap_x1","xmap_x1c","xmap_x2","xmap_x2c","xmap_x3","xmap_x3c"],\
                                    "calc":["arr_absorp1","arr_absorp2","arr_absorp3"]})
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
                                    {"I0":"zap_iodet",\
                                    "It":"zap_idet",\
                                    "XRF":["xmap_x1","xmap_x1c","xmap_x2","xmap_x2c","xmap_x3","xmap_x3c","xmap_icr","xmap_ocr"]})
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
        info["counterdict"] = info.get("counterdict",\
                                    {"I0":"zap_p201_IC",\
                                    "It":"zap_p201_It"})
        info["metadata"] = info.get("metadata","xia")
        super(ESRF_ID16B,self).__init__(**info)

        
factory = InstrumentInfo.factory


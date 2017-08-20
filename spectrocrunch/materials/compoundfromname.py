# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
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

from .compound import compound
from .compoundfromlist import compoundfromlist
from .compoundfromformula import compoundfromformula
from .mixture import mixture
from .types import fractionType
import xraylib

compounddb = {}
compounddb["vacuum"] = compoundfromlist([],[],fractionType.mole,0,name="vacuum")

# triglycerides:
compounddb["trilinolein"] = compoundfromformula("C57H98O6",0.925,name="trilinolein")
compounddb["trilinolenin"] = compoundfromformula("C57H92O6",0.946,name="trilinolenin")
compounddb["triolein"] = compoundfromformula("C57H104O6",0.95,name="triolein")
compounddb["tripalmitin"] = compoundfromformula("C51H98O6",0.8752,name="tripalmitin")
compounddb["tristearin"] = compoundfromformula("C57H110O6",0.909,name="tristearin")
compounddb["tripalmitolein"] = compoundfromformula("C51H92O6",0.929,name="tripalmitolein")
compounddb["triarachidin"] = compoundfromformula("C63H122O6",0.9540,name="triarachidin")
#compounddb["eicosenoin"] = compoundfromformula("C23H44O4",0.9540,name="eicosenoin")

# oils:
# http://www.journal-of-agroalimentary.ro/admin/articole/61602L07_Popa_Vol.18%282%29_2012.pdf
# linolenic (53.21%), oleic (18.51%), linoleic (17.25%), palmitic (6.58 %) and stearic (4.43%).
# https://link.springer.com/content/pdf/10.1007%2Fs11746-004-0879-6.pdf
# linolenic (58.5%), oleic (16.1%), linoleic (14.7%), palmitic (7 %) and stearic (2.9%).
compounddb["linseedoil"] = mixture([compounddb["trilinolenin"],\
                                    compounddb["triolein"],\
                                    compounddb["trilinolein"],\
                                    compounddb["tripalmitin"],\
                                    compounddb["tristearin"]],\
                                    #[0.5321,0.1851,0.1725,0.0658,0.00443],\
                                    [0.585,0.161,0.147,0.07,0.029],\
                                    fractionType.mole).tocompound("linseed oil")

# organics:
compounddb["cellulose"] = compoundfromformula("C6H10O5",1.5,name="")

# 
compounddb["diamond"] = compoundfromformula("C",3.51,name="")

compounddb["silica"] = compoundfromformula("SiO2",2.648,name="silica")
compounddb["quartz"] = compoundfromformula("SiO2",2.66,name="quartz")

compounddb["calcite"] = compoundfromformula("CaCO3",2.7102,name="calcite")
compounddb["gypsum"] = compoundfromformula("CaSO6H4",2.31,name="gypsum")
compounddb["lazurite"] = compoundfromformula("Na3CaAl3Si3O12S",2.4,name="lazurite")
compounddb["hydroxyapatite"] = compoundfromformula("Ca5(PO4)3(OH)",3.,name="hydroxyapatite")

compounddb["rutile"] = compoundfromformula("TiO2",4.25,name="rutile")

compounddb["cobaltblue"] = compoundfromformula("CoAl2O4",8.9,name="cobaltblue")

compounddb["hematite"] = compoundfromformula("Fe2O3",5.26,name="hematite")
compounddb["magnetite"] = compoundfromformula("Fe2O3",5.15,name="magnetite")
compounddb["goethite"] = compoundfromformula("FeOOH",3.8,name="magnetite")
compounddb["vivianite"] = compoundfromformula("Fe3P2O12H16",2.65,name="")
compounddb["prussianblue"] = compoundfromformula("C18Fe7N18",1.83,name="prussian blue")
compounddb["potassiumferricyanide"] = compoundfromformula("K3Fe(CN)6",1.89,name="potassium ferricyanide")
compounddb["potassiumferrocyanide"] = compoundfromformula("K4Fe(CN)6(H2O)3",1.85,name="potassium ferrocyanide")
compounddb["ferrousoxalate"] = compoundfromformula("FeC2O4(H2O)2",2.28,name="ferrous oxalate")
compounddb["ferricoxalate"] = compoundfromformula("Fe2(C2O4)3(H2O)6",2.28,name="ferric oxalate")
compounddb["ferroussulfate"] = compoundfromformula("FeSO4(H2O)7",1.895,name="ferrous sulfate")
compounddb["ferricsulfate"] = compoundfromformula("Fe2(SO4)3(H2O)5",1.898,name="ferric sulfate")
compounddb["ferricchloride"] = compoundfromformula("FeCl3(H2O)6",1.82,name="ferric chloride")

compounddb["cadmiumsulfide"] = compoundfromformula("CdS",4.826,name="cadmiumsulfide")

compounddb["corderoite"] = compoundfromformula("Hg3S2Cl2",6.845,name="corderoite")

compounddb["hydrocerussite"] = compoundfromformula("Pb3C2O8H2",6.8,name="hydrocerussite")
compounddb["leadchromate"] = compoundfromformula("PbCrO4",6.3,name="leadchromate")

# tapes/foils:
data = xraylib.GetCompoundDataNISTByName("Kapton Polyimide Film")
compounddb["kapton"] = compound(data["Elements"],data["massFractions"],fractionType.weight,data["density"],name="kapton")
compounddb["ultralene"] = compound(data["Elements"],data["massFractions"],fractionType.weight,data["density"],name="ultralene")

# resins/plastics:
compounddb["pmma"] = compoundfromformula("C5O2H8",1.18,name="polymethyl methacrylate ")
compounddb["pp"] = compoundfromformula("C3H6",0.86,name="polypropylene")
compounddb["pe"] = compoundfromformula("C2H4",0.95,name="polyethylene")
compounddb["pet"] = compoundfromformula("C10H8O4",1.38,name="polyethylene terephthalate ")
compounddb["pan"] = compoundfromformula("C3H3N",1.184,name="polyacrylonitrile")

def compoundfromname(name):
    return compounddb[name]





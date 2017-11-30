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

from . import compound
from . import compoundfromlist
from . import compoundfromformula
from . import mixture
from .types import fractionType
import xraylib

compounddb = {}
compounddb["vacuum"] = compoundfromlist.CompoundFromList([],[],fractionType.mole,0,name="vacuum")

# triglycerides:
compounddb["trilinolein"] = compoundfromformula.CompoundFromFormula("C57H98O6",0.925,name="trilinolein")
compounddb["trilinolenin"] = compoundfromformula.CompoundFromFormula("C57H92O6",0.946,name="trilinolenin")
compounddb["triolein"] = compoundfromformula.CompoundFromFormula("C57H104O6",0.95,name="triolein")
compounddb["tripalmitin"] = compoundfromformula.CompoundFromFormula("C51H98O6",0.8752,name="tripalmitin")
compounddb["tristearin"] = compoundfromformula.CompoundFromFormula("C57H110O6",0.909,name="tristearin")
compounddb["tripalmitolein"] = compoundfromformula.CompoundFromFormula("C51H92O6",0.929,name="tripalmitolein")
compounddb["triarachidin"] = compoundfromformula.CompoundFromFormula("C63H122O6",0.9540,name="triarachidin")
#compounddb["eicosenoin"] = compoundfromformula.CompoundFromFormula("C23H44O4",0.9540,name="eicosenoin")

# oils:
# http://www.journal-of-agroalimentary.ro/admin/articole/61602L07_Popa_Vol.18%282%29_2012.pdf
# linolenic (53.21%), oleic (18.51%), linoleic (17.25%), palmitic (6.58 %) and stearic (4.43%).
# https://link.springer.com/content/pdf/10.1007%2Fs11746-004-0879-6.pdf
# linolenic (58.5%), oleic (16.1%), linoleic (14.7%), palmitic (7 %) and stearic (2.9%).
compounddb["linseed oil"] = mixture.Mixture([compounddb["trilinolenin"],\
                                    compounddb["triolein"],\
                                    compounddb["trilinolein"],\
                                    compounddb["tripalmitin"],\
                                    compounddb["tristearin"]],\
                                    #[0.5321,0.1851,0.1725,0.0658,0.00443],\
                                    [0.585,0.161,0.147,0.07,0.029],\
                                    fractionType.mole).tocompound("linseed oil")

# organics:
compounddb["cellulose"] = compoundfromformula.CompoundFromFormula("C6H10O5",1.5,name="")

# low-Z minerals
compounddb["diamond"] = compoundfromformula.CompoundFromFormula("C",3.51,name="")
compounddb["silica"] = compoundfromformula.CompoundFromFormula("SiO2",2.648,name="silica")
compounddb["quartz"] = compoundfromformula.CompoundFromFormula("SiO2",2.66,name="quartz")

# Ca minerals
compounddb["calcite"] = compoundfromformula.CompoundFromFormula("CaCO3",2.7102,name="calcite")
compounddb["gypsum"] = compoundfromformula.CompoundFromFormula("CaSO6H4",2.31,name="gypsum")
compounddb["lazurite"] = compoundfromformula.CompoundFromFormula("Na3CaAl3Si3O12S",2.4,name="lazurite")
compounddb["hydroxyapatite"] = compoundfromformula.CompoundFromFormula("Ca5(PO4)3(OH)",3.,name="hydroxyapatite")

# Ti minerals
compounddb["rutile"] = compoundfromformula.CompoundFromFormula("TiO2",4.25,name="rutile")

# Co pigments
compounddb["cobalt blue"] = compoundfromformula.CompoundFromFormula("CoAl2O4",8.9,name="cobalt blue")

# Fe minerals
compounddb["hematite"] = compoundfromformula.CompoundFromFormula("Fe2O3",5.26,name="hematite")
compounddb["magnetite"] = compoundfromformula.CompoundFromFormula("Fe2O3",5.15,name="magnetite")
compounddb["goethite"] = compoundfromformula.CompoundFromFormula("FeOOH",3.8,name="magnetite")
compounddb["vivianite"] = compoundfromformula.CompoundFromFormula("Fe3P2O12H16",2.65,name="")
compounddb["prussian blue"] = compoundfromformula.CompoundFromFormula("C18Fe7N18",1.83,name="prussian blue")
compounddb["potassium ferricyanide"] = compoundfromformula.CompoundFromFormula("K3Fe(CN)6",1.89,name="potassium ferricyanide")
compounddb["potassium ferrocyanide"] = compoundfromformula.CompoundFromFormula("K4Fe(CN)6(H2O)3",1.85,name="potassium ferrocyanide")
compounddb["ferrous oxalate"] = compoundfromformula.CompoundFromFormula("FeC2O4(H2O)2",2.28,name="ferrous oxalate")
compounddb["ferric oxalate"] = compoundfromformula.CompoundFromFormula("Fe2(C2O4)3(H2O)6",2.28,name="ferric oxalate")
compounddb["ferrous sulfate"] = compoundfromformula.CompoundFromFormula("FeSO4(H2O)7",1.895,name="ferrous sulfate")
compounddb["ferric sulfate"] = compoundfromformula.CompoundFromFormula("Fe2(SO4)3(H2O)5",1.898,name="ferric sulfate")
compounddb["ferric chloride"] = compoundfromformula.CompoundFromFormula("FeCl3(H2O)6",1.82,name="ferric chloride")

# Cd pigments
compounddb["cadmium sulfide"] = compoundfromformula.CompoundFromFormula("CdS",4.826,name="cadmium sulfide")

# Hg pigments
compounddb["corderoite"] = compoundfromformula.CompoundFromFormula("Hg3S2Cl2",6.845,name="corderoite")

# Pb pigments
compounddb["hydrocerussite"] = compoundfromformula.CompoundFromFormula("Pb3C2O8H2",6.8,name="hydrocerussite")
compounddb["lead chromate"] = compoundfromformula.CompoundFromFormula("PbCrO4",6.3,name="lead chromate")

# Cu pigments
compounddb["copper acetate"] = compoundfromformula.CompoundFromFormula("Cu(CH3CO2)2", 1.882,name="copper acetate")
compounddb["verdigris"] = compoundfromformula.CompoundFromFormula("Cu2(CH3CO2)4(H2O)2", 1.882,name="verdigris")
compounddb["paris green"] = compoundfromformula.CompoundFromFormula("Cu2(CH3CO2)4(H2O)2", 1.1,name="paris green")

# tapes/foils:
data = xraylib.GetCompoundDataNISTByName("Kapton Polyimide Film")
compounddb["kapton"] = compound.Compound(data["Elements"],data["massFractions"],fractionType.weight,data["density"],name="kapton")
#compounddb["ultralene"] = compound.Compound(data["Elements"],data["massFractions"],fractionType.weight,data["density"],name="ultralene")
compounddb["mylar"] = compoundfromformula.CompoundFromFormula("C10H8O4",1.38,name="mylar")
compounddb["ultralene"] = compoundfromformula.CompoundFromFormula("C10H8O4",1.38,name="ultralene")

# polymers:
compounddb["pmma"] = compoundfromformula.CompoundFromFormula("C5O2H8",1.18,name="pmma")
compounddb["pp"] = compoundfromformula.CompoundFromFormula("C3H6",0.86,name="pp")
compounddb["pe"] = compoundfromformula.CompoundFromFormula("C2H4",0.95,name="pe")
compounddb["pet"] = compoundfromformula.CompoundFromFormula("C10H8O4",1.38,name="pet")
compounddb["pan"] = compoundfromformula.CompoundFromFormula("C3H3N",1.184,name="pan")
compounddb["pva"] = compoundfromformula.CompoundFromFormula("C3H3N",1.19,name="pva")

# tape (adhesive on a plastic)
compounddb["sulfur-free tape"] = mixture.Mixture([compounddb["pva"],compounddb["pe"]],\
                                    [0.5,0.5],fractionType.mole).tocompound("sulfur-free tape")

# windows
compounddb["silicon nitride"] = compoundfromformula.CompoundFromFormula("Si3N4",3.44,name="silicon nitride")

def compoundfromname(name):
    return compounddb[name]





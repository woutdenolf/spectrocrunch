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

import warnings
from . import compound
from . import compoundfromlist
from . import compoundfromformula
from . import mixture
from . import types
from ..utils import instance

compounddb = {}

try:
    import xraylib
except ImportError:
    xraylib = None
    warnings.warn("xraylib is not installed", ImportWarning)
else:
    compounddb["vacuum"] = compoundfromlist.CompoundFromList(
        [], [], types.fraction.mole, 0, name="vacuum")
    compounddb["dummy"] = compoundfromlist.CompoundFromList(
        ['H'], [1], types.fraction.mole, 0, name="dummy")

    # triglycerides:
    compounddb["trilinolein"] = compoundfromformula.CompoundFromFormula(
        "C57H98O6", 0.925, name="trilinolein")
    compounddb["trilinolenin"] = compoundfromformula.CompoundFromFormula(
        "C57H92O6", 0.946, name="trilinolenin")
    compounddb["triolein"] = compoundfromformula.CompoundFromFormula(
        "C57H104O6", 0.95, name="triolein")
    compounddb["tripalmitin"] = compoundfromformula.CompoundFromFormula(
        "C51H98O6", 0.8752, name="tripalmitin")
    compounddb["tristearin"] = compoundfromformula.CompoundFromFormula(
        "C57H110O6", 0.909, name="tristearin")
    compounddb["tripalmitolein"] = compoundfromformula.CompoundFromFormula(
        "C51H92O6", 0.929, name="tripalmitolein")
    compounddb["triarachidin"] = compoundfromformula.CompoundFromFormula(
        "C63H122O6", 0.9540, name="triarachidin")
    #compounddb["eicosenoin"] = compoundfromformula.CompoundFromFormula("C23H44O4",0.9540,name="eicosenoin")

    # oils:
    # http://www.journal-of-agroalimentary.ro/admin/articole/61602L07_Popa_Vol.18%282%29_2012.pdf
    # linolenic (53.21%), oleic (18.51%), linoleic (17.25%), palmitic (6.58 %) and stearic (4.43%).
    # https://link.springer.com/content/pdf/10.1007%2Fs11746-004-0879-6.pdf
    # linolenic (58.5%), oleic (16.1%), linoleic (14.7%), palmitic (7 %) and stearic (2.9%).
    compounddb["linseed oil"] = mixture.Mixture([compounddb["trilinolenin"],
                                                 compounddb["triolein"],
                                                 compounddb["trilinolein"],
                                                 compounddb["tripalmitin"],
                                                 compounddb["tristearin"]],\
                                                # [0.5321,0.1851,0.1725,0.0658,0.00443],\
                                                [0.585, 0.161, 0.147, 0.07, 0.029],\
                                                types.fraction.mole).tocompound("linseed oil")

    # organics:
    compounddb["cellulose"] = compoundfromformula.CompoundFromFormula(
        "C6H10O5", 1.5, name="cellulose")
    # http://dx.doi.org/10.1100/tsw.2009.27:
    compounddb["human hair"] = compoundfromlist.CompoundFromList(['C', 'O', 'N', 'H', 'S'],
                                                                 [0.49, 0.30, 0.145,
                                                                     0.035, 0.03],
                                                                 types.fraction.mass, 1.33, name="human hair")
    compounddb["hair"] = compounddb["human hair"]

    # low-Z minerals
    compounddb["diamond"] = compoundfromformula.CompoundFromFormula(
        "C", 3.51, name="")
    compounddb["Na2O"] = compoundfromformula.CompoundFromFormula(
        "Na2O", 2.27, name="Na2O")
    compounddb["silica"] = compoundfromformula.CompoundFromFormula(
        "SiO2", 2.169, name="silica")
    compounddb["quartz"] = compoundfromformula.CompoundFromFormula(
        "SiO2", 2.648, name="quartz")
    compounddb["corundum"] = compoundfromformula.CompoundFromFormula(
        "Al2O3", 4.02, name="corundum")
    compounddb["alumina"] = compoundfromformula.CompoundFromFormula(
        "Al2O3", 3.987, name="alumina")
    compounddb["periclase"] = compoundfromformula.CompoundFromFormula(
        "MgO", 3.6, name="periclase")
    compounddb["K2O"] = compoundfromformula.CompoundFromFormula(
        "K2O", 2.2, name="K2O")
    compounddb["P2O5"] = compoundfromformula.CompoundFromFormula(
        "P2O5", 2.39, name="P2O5")
    compounddb["NaCl"] = compoundfromformula.CompoundFromFormula(
        "NaCl", 2.16, name="NaCl")
    compounddb["KCl"] = compoundfromformula.CompoundFromFormula(
        "KCl", 1.984, name="KCl")

    # Ca minerals
    compounddb["calcite"] = compoundfromformula.CompoundFromFormula(
        "CaCO3", 2.7102, name="calcite")
    compounddb["gypsum"] = compoundfromformula.CompoundFromFormula(
        "CaSO6H4", 2.31, name="gypsum")
    compounddb["lazurite"] = compoundfromformula.CompoundFromFormula(
        "Na3CaAl3Si3O12S", 2.4, name="lazurite")
    compounddb["hydroxyapatite"] = compoundfromformula.CompoundFromFormula(
        "Ca5(PO4)3(OH)", 3., name="hydroxyapatite")
    compounddb["lime"] = compoundfromformula.CompoundFromFormula(
        "CaO", 3.34, name="lime")

    # Ti minerals
    compounddb["rutile"] = compoundfromformula.CompoundFromFormula(
        "TiO2", 4.25, name="rutile")
    compounddb["anatase"] = compoundfromformula.CompoundFromFormula(
        "TiO2", 4.23, name="anatase")

    # Co pigments
    compounddb["cobalt blue"] = compoundfromformula.CompoundFromFormula(
        "CoAl2O4", 8.9, name="cobalt blue")

    # Fe minerals
    compounddb["hematite"] = compoundfromformula.CompoundFromFormula(
        "Fe2O3", 5.26, name="hematite")
    compounddb["magnetite"] = compoundfromformula.CompoundFromFormula(
        "Fe2O3", 5.15, name="magnetite")
    compounddb["ferric oxide"] = compoundfromformula.CompoundFromFormula(
        "Fe2O3", 5.245, name="ferric oxide")
    compounddb["ferrous oxide"] = compoundfromformula.CompoundFromFormula(
        "FeO", 5.745, name="ferrous oxide")
    compounddb["goethite"] = compoundfromformula.CompoundFromFormula(
        "FeOOH", 3.8, name="goethite")
    compounddb["vivianite"] = compoundfromformula.CompoundFromFormula(
        "Fe3(PO4)2(H2O)8", 2.65, name="vivianite")
    compounddb["prussian blue"] = compoundfromformula.CompoundFromFormula(
        "C18Fe7N18", 1.83, name="prussian blue")
    compounddb["potassium ferricyanide"] = compoundfromformula.CompoundFromFormula(
        "K3Fe(CN)6", 1.89, name="potassium ferricyanide")
    compounddb["potassium ferrocyanide"] = compoundfromformula.CompoundFromFormula(
        "K4Fe(CN)6(H2O)3", 1.85, name="potassium ferrocyanide")
    compounddb["ferrous oxalate"] = compoundfromformula.CompoundFromFormula(
        "FeC2O4(H2O)2", 2.28, name="ferrous oxalate")
    compounddb["ferric oxalate"] = compoundfromformula.CompoundFromFormula(
        "Fe2(C2O4)3(H2O)6", 2.28, name="ferric oxalate")
    compounddb["ferrous sulfate"] = compoundfromformula.CompoundFromFormula(
        "FeSO4(H2O)7", 1.895, name="ferrous sulfate")
    compounddb["ferric sulfate"] = compoundfromformula.CompoundFromFormula(
        "Fe2(SO4)3(H2O)5", 1.898, name="ferric sulfate")
    compounddb["ferric chloride"] = compoundfromformula.CompoundFromFormula(
        "FeCl3(H2O)6", 1.82, name="ferric chloride")

    # Cr minerals
    compounddb["chromia"] = compoundfromformula.CompoundFromFormula(
        "Cr2O3", 5.22, name="chromia")

    # Ni minerals
    compounddb["bunsenite"] = compoundfromformula.CompoundFromFormula(
        "NiO", 6.67, name="bunsenite")

    # Mn minerals
    compounddb["MnO"] = compoundfromformula.CompoundFromFormula(
        "MnO", 5.745, name="MnO")

    # Cd pigments
    compounddb["cadmium sulfide"] = compoundfromformula.CompoundFromFormula(
        "CdS", 4.826, name="cadmium sulfide")

    # Hg pigments
    compounddb["corderoite"] = compoundfromformula.CompoundFromFormula(
        "Hg3S2Cl2", 6.845, name="corderoite")

    # Pb pigments
    compounddb["cerussite"] = compoundfromformula.CompoundFromFormula(
        "PbCO3", 6.58, name="cerussite")
    compounddb["hydrocerussite"] = compoundfromformula.CompoundFromFormula(
        "Pb3C2O8H2", 6.8, name="hydrocerussite")
    compounddb["lead chromate"] = compoundfromformula.CompoundFromFormula(
        "PbCrO4", 6.3, name="lead chromate")

    # Cu pigments
    compounddb["copper acetate"] = compoundfromformula.CompoundFromFormula(
        "Cu(CH3CO2)2", 1.882, name="copper acetate")
    compounddb["verdigris"] = compoundfromformula.CompoundFromFormula(
        "Cu2(CH3CO2)4(H2O)2", 1.882, name="verdigris")
    compounddb["paris green"] = compoundfromformula.CompoundFromFormula(
        "Cu2(CH3CO2)4(H2O)2", 1.1, name="paris green")
    compounddb["malachite"] = compoundfromformula.CompoundFromFormula(
        "Cu2CO3(OH)2", 4.03, name="malachite")
    compounddb["cuprorivaite"] = compoundfromformula.CompoundFromFormula(
        "CaCuSi4O10", 3.1, name="cuprorivaite")
    compounddb["tenorite"] = compoundfromformula.CompoundFromFormula(
        "CuO", 6.45, name="tenorite")
    compounddb["dioptase"] = compoundfromformula.CompoundFromFormula(
        "CuSiO2(OH)2", 3.3, name="dioptase")
    compounddb["atacamite"] = compoundfromformula.CompoundFromFormula(
        "Cu2Cl(OH)3", 3.756, name="atacamite")

    # polymers:
    compounddb["pmma"] = compoundfromformula.CompoundFromFormula(
        "C5O2H8", 1.18, name="pmma")
    compounddb["pe"] = compoundfromformula.CompoundFromFormula(
        "C2H4", 0.95, name="pe")
    compounddb["pan"] = compoundfromformula.CompoundFromFormula(
        "C3H3N", 1.184, name="pan")
    compounddb["pva"] = compoundfromformula.CompoundFromFormula(
        "C3H3N", 1.19, name="pva")

    compounddb["pp"] = compoundfromformula.CompoundFromFormula(
        "C3H6", 0.86, name="pp")
    compounddb["specx pp"] = compoundfromformula.CompoundFromFormula(
        "C3H6", 1.1289633445, name="pp")

    data = xraylib.GetCompoundDataNISTByName("Kapton Polyimide Film")
    compounddb["kapton"] = compound.Compound(
        data["Elements"], data["massFractions"], types.fraction.mass, data["density"], name="kapton")
    compounddb["specx kapton"] = compound.Compound(
        data["Elements"], data["massFractions"], types.fraction.mass, 1.35297609903, name="kapton")

    compounddb["pet"] = compoundfromformula.CompoundFromFormula(
        "C10H8O4", 1.38, name="pet")
    compounddb["mylar"] = compoundfromformula.CompoundFromFormula(
        "C10H8O4", 1.38, name="mylar")  # same as pet
    compounddb["specx mylar"] = compoundfromformula.CompoundFromFormula(
        "C10H8O4", 0.933152791636, name="mylar")

    compounddb["pc"] = compoundfromformula.CompoundFromFormula(
        "C15H16O2", 1.2, name="pc")
    compounddb["specx pc"] = compoundfromformula.CompoundFromFormula(
        "C15H16O2", 0.485864569962, name="pc")
    compounddb["specx ultralene"] = compoundfromformula.CompoundFromFormula(
        "C15H16O2", 0.452164634083, name="ultralene")
    compounddb["ultralene"] = compounddb["specx ultralene"]  # 4.064e-4 cm

    compounddb["moxtek ap3.3"] = compoundfromlist.CompoundFromList(['B', 'C', 'N', 'O', 'Al'],
                                                                   [0.13336917388076333, 0.5117702789499711, 0.11331306454114157,
                                                                       0.20186742974474034, 0.039680052883383735],
                                                                   types.fraction.mole, 0.757636153465, name="moxtek ap3.3")
    compounddb["moxtek ap3.7"] = compoundfromlist.CompoundFromList(['B', 'C', 'N', 'O', 'Al'],
                                                                   [0.14109395570244615, 0.3969446446586894, 0.20446566594853216,
                                                                       0.195183216894773, 0.06231251679555934],
                                                                   types.fraction.mole, 1.23845839755, name="moxtek ap3.7")

    compounddb["epoxy resin"] = compoundfromformula.CompoundFromFormula(
        "C21H25ClO5", 2., name="epoxy resin")

    # tape (adhesive on plastic)
    compounddb["sulfur-free tape"] = mixture.Mixture([compounddb["pva"], compounddb["pe"]],
                                                     [0.5, 0.5], types.fraction.mole).tocompound("sulfur-free tape")  # 50 um

    # windows
    compounddb["silicon nitride"] = compoundfromformula.CompoundFromFormula(
        "Si3N4", 3.44, name="silicon nitride")

    # rocks
    compounddb["granite"] = mixture.Mixture([compounddb["silica"],
                                             compounddb["alumina"],
                                             compounddb["K2O"],
                                             compounddb["Na2O"],
                                             compounddb["lime"],
                                             compounddb["ferrous oxide"],
                                             compounddb["ferric oxide"],
                                             compounddb["periclase"],
                                             compounddb["rutile"],
                                             compounddb["P2O5"],
                                             compounddb["MnO"]],
                                            [0.7204, 0.1442, 0.0412, 0.0369, 0.0182, 0.0168,
                                                0.0122, 0.0071, 0.0030, 0.0012, 0.0005],
                                            types.fraction.mole).tocompound("granite")

    compounddb["orthoclase"] = compoundfromformula.CompoundFromFormula(
        "KAlSi3O8", 2.56, name="orthoclase")
    compounddb["albite"] = compoundfromformula.CompoundFromFormula(
        "NaAlSi3O8", 2.62, name="albite")
    compounddb["anorthite"] = compoundfromformula.CompoundFromFormula(
        "CaAl2Si2O8", 2.73, name="anorthite")

    # other
    compounddb["ice"] = compoundfromformula.CompoundFromFormula(
        "H2O", 0.9167, name="ice")
    compounddb["water"] = compoundfromformula.CompoundFromFormula(
        "H2O", 0.9998, name="water")

    compounddb["ecoli dry"] = compoundfromlist.CompoundFromList(['C', 'O', 'N', 'H', 'P', 'S', 'K', 'Na', 'Ca', 'Mg', 'Cl', 'Fe'],
                                                                [50, 20, 14, 8, 3, 1, 1, 1, 0.5, 0.5, 0.5, 0.5], types.fraction.mass, 1.3, name="ecoli dry")

    compounddb["ecoli"] = mixture.Mixture([compounddb["water"], compounddb["ecoli dry"]],
                                          [0.7, 0.3], types.fraction.mole).tocompound("ecoli")

    data = xraylib.GetCompoundDataNISTByName("Air, Dry (near sea level)")
    compounddb["air dry"] = compound.Compound(
        data["Elements"], data["massFractions"], types.fraction.mass, data["density"], name="air dry")
    compounddb["air"] = compounddb["air dry"]

registry = compounddb.keys()


def factory(name):
    return compounddb[name]


def search(name):
    name = name.lower()
    ret = [k for k in registry if name in k.lower()]
    if len(ret) > 1:
        ret2 = [k for k in registry if name == k.lower()]
        if ret2:
            ret = ret2
    if ret:
        return ret[0]
    else:
        return None


def compoundfromname(name):
    return factory(name)

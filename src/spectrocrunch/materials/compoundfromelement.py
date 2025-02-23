import re
from . import compoundfromformula
from . import element

db = {}
db["oxide"] = {
    "H": {1: ("H2O", 0.997)},
    "Li": {1: ("Li2O", 2.01)},
    "Be": {2: ("BeO", 3.02)},
    "B": {3: ("B2O3", 2.46)},
    "C": {4: ("CO2", 1.98)},
    "Na": {1: ("Na2O", 2.27)},
    "Mg": {2: ("MgO", 3.58)},
    "Al": {3: ("Al2O3", 3.95)},
    "Si": {4: ("SiO2", 2.65)},
    "P": {5: ("P2O5", 2.39)},
    "S": {4: ("SO2", 2.63), 6: ("SO3", 1.92)},
    "K": {1: ("K2O", 2.35)},
    "Ca": {2: ("CaO", 3.35)},
    "Sc": {3: ("Sc2O3", 3.86)},
    "Ti": {3: ("Ti2O3", 4.49), 4: ("TiO2", 4.23)},
    "V": {5: ("V2O5", 3.36)},
    "Cr": {3: ("Cr2O3", 5.22)},
    "Mn": {2: ("MnO", 5.37), 8 / 3.0: ("Mn3O4", 4.86), 4: ("MnO2", 5.03)},
    "Fe": {2: ("FeO", 5.74), 3: ("Fe2O3", 5.24)},
    "Co": {2: ("CoO", 6.44)},
    "Ni": {2: ("NiO", 6.67)},
    "Cu": {2: ("CuO", 6.31)},
    "Zn": {2: ("ZnO", 5.61)},
    "Ga": {3: ("Ga2O3", 6.44)},
    "Ge": {4: ("GeO2", 4.25)},
    "As": {3: ("As2O3", 3.74)},
    "Rb": {1: ("Rb2O", 4.0)},
    "Sr": {2: ("SrO", 4.7)},
    "Y": {3: ("Y2O3", 5.01)},
    "Zr": {4: ("ZrO2", 5.68)},
    "Nb": {5: ("Nb2O5", 4.6)},
    "Mo": {6: ("MoO3", 4.69)},
    "Sn": {4: ("SnO2", 6.95)},
    "Sb": {3: ("Sb2O3", 5.2)},
    "Cs": {1: ("Cs2O", 4.65)},
    "Ba": {2: ("BaO", 5.72)},
    "La": {3: ("La2O3", 6.51)},
    "Ce": {3: ("Ce2O3", 6.2), 4: ("CeO2", 7.22)},
    "Pr": {3: ("Pr2O3", 6.9)},
    "Nd": {3: ("Nd2O3", 7.24)},
    "Sm": {3: ("Sm2O3", 8.35)},
    "Eu": {2: ("EuO", 8.2), 3: ("Eu2O3", 7.4)},
    "Gd": {3: ("Gd2O3", 7.41)},
    "Tb": {3: ("Tb2O3", 7.9)},
    "Dy": {3: ("Dy2O3", 7.8)},
    "Ho": {3: ("Ho2O3", 8.41)},
    "Er": {3: ("Er2O3", 8.64)},
    "Tm": {3: ("Tm2O3", 8.6)},
    "Yb": {3: ("Yb2O3", 9.17)},
    "Lu": {3: ("Lu2O3", 9.42)},
    "Hf": {4: ("HfO2", 9.68)},
    "Ta": {5: ("Ta2O5", 8.2)},
    "W": {6: ("WO3", 7.16)},
    "Au": {3: ("Au2O3", 11.34)},
    "Pb": {2: ("PbO", 9.53)},
    "Th": {4: ("ThO2", 10.0)},
    "U": {4: ("UO2", 10.97), 16 / 3.0: ("U3O8", 8.39)},
}

db["sulfide"] = {
    "Fe": {1: ("Fe2S", 4.7), 2: ("FeS", 4.84)},
    "Cu": {1: ("Cu2S", 5.6), 2: ("CuS", 4.76)},
    "Ni": {2: ("NiS", 5.87)},
    "Zn": {2: ("ZnS", 4.09)},
    "Ag": {1: ("Ag2S", 7.23)},
    "Pb": {2: ("PbS", 7.6)},
    "Hg": {2: ("HgS", 8.1)},
    "As": {2: ("AsS", 3.56), 3: ("As2S3", 3.43)},
    "Sb": {3: ("Sb2S3", 6.5)},
    "Mo": {4: ("MoS2", 5.06)},
    "Cd": {2: ("CdS", 4.826)},
}

# TODO:
# https://en.wikipedia.org/wiki/Sulfate
# http://www.endmemo.com/chem/common/chloride.php


class CompoundFromElement(compoundfromformula.CompoundFromFormula):
    def __init__(
        self, el, oxidationstate=None, typ="oxide", defaultoxidationstate="max"
    ):
        if typ not in db:
            raise RuntimeError(
                "Type {} does not exist (valid types: {})".format(typ, db.keys())
            )
        if el not in db[typ]:
            raise RuntimeError(
                "A compound with element {} of type {} does not exist".format(el, typ)
            )
        cmpds = db[typ][el]

        if oxidationstate is None:
            if len(cmpds) == 1:
                cmpd = cmpds.values()[0]
            else:
                if defaultoxidationstate == "max":
                    func = max
                elif defaultoxidationstate == "min":
                    func = min
                else:
                    raise RuntimeError(
                        "Select one of the available {} oxidations states: {}".format(
                            el, cmpds.keys()
                        )
                    )
                oxidationstate = func(cmpds.keys())
                cmpd = cmpds[oxidationstate]
        elif oxidationstate in cmpds:
            cmpd = cmpds[oxidationstate]
        else:
            raise RuntimeError(
                "Oxidation state {} does not exist for {} {}".format(
                    oxidationstate, el, typ
                )
            )

        formula, density = cmpd
        super(CompoundFromElement, self).__init__(formula, density=density)

        self.oxidationstate = oxidationstate

        parse = re.compile(
            r"^(?P<Z>[A-Z][a-z]?)(?P<a>\d?)(?P<X>[A-Z][a-z]?)(?P<b>\d?)$"
        )
        self.formula = parse.match(formula).groupdict()

        self.formula["Z"] = element.Element(self.formula["Z"])
        self.formula["X"] = element.Element(self.formula["X"])
        for k in ["a", "b"]:
            if self.formula[k]:
                self.formula[k] = float(self.formula[k])
            else:
                self.formula[k] = 1.0

    def _massfraction_conv(self):
        # Element Z is present only as C = ZaXb (X is typically oxigen)
        # wi = ni.MMi / sum_j(nj.MMj)
        # wCi/wZi = nCi.MMCi/(nZi.MMZi)
        # wCi = MMCi/MMZi . nCi/nZi . wZi
        #     = MMCi/(aZi.MMZi)  .wZi
        MMC = self.molarmass()
        MMZ = element.Element(self.formula["Z"]).molarmass()
        aZ = float(self.formula["a"])
        return MMC / (aZ * MMZ)

    def massfraction_element_to_compound(self, wfrac):
        return wfrac * self._massfraction_conv()

    def massfraction_compound_to_element(self, wfrac):
        return wfrac / self._massfraction_conv()


class Oxide(CompoundFromElement):
    def __init__(self, *args, **kwargs):
        kwargs["typ"] = "oxide"
        super(Oxide, self).__init__(*args, **kwargs)


class Sulfide(CompoundFromElement):
    def __init__(self, *args, **kwargs):
        kwargs["typ"] = "sulfide"
        super(Oxide, self).__init__(*args, **kwargs)

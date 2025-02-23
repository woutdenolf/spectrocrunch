from . import multielementbase
from . import element
from . import types
from . import stoichiometry
from ..utils import instance

import numpy as np
import fisx


class Compound(multielementbase.MultiElementBase):
    """Interface to a compound"""

    def __init__(
        self, elements, frac, fractype=None, density=None, nrefrac=1, name=None
    ):
        """
        Args:
            elements(list): list of elements (["Fe","O"] or [element("Fe"),element("O")])
            frac(list[float]): element fractions
            fractype(types.fraction): element fraction type
            density(num): compound density in g/cm^3
            nrefrac(num): refractive index
            name(Optional[str]): compound name
        """

        # Compound name
        if name is None:
            self.name = str(id(self))
        else:
            self.name = name
        self.nrefrac = float(nrefrac)

        # Element mole fractions
        if fractype == types.fraction.mole:
            nfrac = frac  # keep unnormalized!
        elif fractype == types.fraction.volume:
            # would be possible if you give the element densities in the compound
            # (which is not the same as the pure element density) but that's an unlikely given
            raise ValueError("Cannot create a compound from elemental volume fractions")
        else:
            elements = [element.Element(e) for e in elements]
            MM = np.asarray([e.MM for e in elements])
            nfrac = stoichiometry.frac_weight_to_mole(
                np.asarray(frac), MM
            )  # normalized

        # Elements (no duplicates)
        self._compose_elements(elements, nfrac)

        # Compound density
        if density == 0 or density is None:
            if len(self._elements) == 1:
                self.density = next(iter(self._elements.keys())).density
            else:
                # rho = [e.density for e in self._elements]
                # self.density = np.mean(rho) # not based on anything, just a value
                if len(self._elements) == 0:
                    self.density = 0.0
                else:
                    self.density = 1.0  # approx. density of water
        else:
            self.density = float(density)

        self.isscatterer = True

    def _compose_elements(self, elements, nfrac):
        self._elements = {}
        for e, n in zip(elements, nfrac):
            if not isinstance(e, element.Element):
                e = element.Element(e)
            if e in self._elements:
                self._elements[e] += float(n)
            else:
                self._elements[e] = float(n)

    def __getstate__(self):
        return {
            "name": self.name,
            "nrefrac": self.nrefrac,
            "density": self.density,
            "elements_keys": list(self._elements.keys()),
            "elements_values": list(self._elements.values()),
            "isscatterer": self.isscatterer,
        }

    def __setstate__(self, state):
        self.name = state["name"]
        self.nrefrac = state["nrefrac"]
        self.density = state["density"]
        self._elements = dict(zip(state["elements_keys"], state["elements_values"]))
        self.isscatterer = state["isscatterer"]

    def addelements(self, elements, fractions, fractype, density=None):
        """Add an element to the compound

        Args:
            elements(str or Element or list): "Fe" or Element("Fe")
            frac(num or list): element fraction
            fractype(types.fraction): element fraction type
            density(Optional(num)): new compound density in g/cm^3
        """
        if fractype == types.fraction.volume:
            raise ValueError("Cannot create a compound from elemental volume fractions")
        if not instance.isarray(elements):
            elements = [elements]
        if not instance.isarray(fractions):
            fractions = [fractions]
        fractions = np.asarray(fractions)
        if fractype == types.fraction.mole:
            nfrac = self.equivalents()
            elements = list(nfrac.keys()) + [element.Element(el) for el in elements]
            nfrac = np.asarray(list(nfrac.values()))
            nfrac = stoichiometry.add_frac(nfrac, fractions)
        else:
            wfrac = self.massfractions()
            elements = list(wfrac.keys()) + [element.Element(el) for el in elements]
            wfrac = np.asarray(list(wfrac.values()))
            wfrac = stoichiometry.add_frac(wfrac, fractions)
            MM = np.asarray([e.MM for e in elements])
            nfrac = stoichiometry.frac_weight_to_mole(wfrac, MM)
        self._compose_elements(elements, nfrac)
        if density:
            self.density = density

    def change_fractions(self, dfrac, fractype):
        """Change the element fractions

        Args:
            dfrac(dict): element fractions
            fractype(types.fraction): fraction type
        """
        elements = list(dfrac.keys())
        fractions = np.asarray(list(dfrac.values()))
        # Element mole fractions
        if fractype == types.fraction.mole:
            nfrac = fractions  # keep unnormalized!
        elif fractype == types.fraction.volume:
            # would be possible if you give the element densities in the compound
            # (which is not the same as the pure element density) but that's an unlikely given
            raise ValueError("Cannot create a compound from elemental volume fractions")
        else:
            MM = np.asarray([e.MM for e in elements])
            nfrac = stoichiometry.frac_weight_to_mole(fractions, MM)  # normalized
        # Elements
        self._compose_elements(elements, nfrac)

    def __getitem__(self, el):
        return self._elements[el]

    @property
    def elements(self):
        return list(self._elements.keys())

    @property
    def parts(self):
        return self._elements

    def molarmass(self):
        nfrac = self.equivalents()
        MM = np.asarray([e.MM for e in nfrac])
        nfrac = np.asarray(list(nfrac.values()))
        return (MM * nfrac).sum()

    def molarmasseff(self):
        nfrac = self.molefractions()
        MM = np.asarray([e.MM for e in nfrac])
        nfrac = np.asarray(list(nfrac.values()))
        return (MM * nfrac).sum()

    @property
    def Zeff(self):
        nfrac = self.molefractions()
        Z = np.asarray([e.Z for e in nfrac])
        nfrac = np.asarray(list(nfrac.values()))
        return (Z * nfrac).sum()

    def molefractions(self):
        nfrac = np.asarray(list(self._elements.values()))
        nfrac /= nfrac.sum()
        return dict(zip(self._elements.keys(), nfrac))

    def equivalents(self):
        return dict(self._elements)

    def massfractions(self):
        nfrac = self.equivalents()
        MM = np.asarray([e.MM for e in nfrac])
        wfrac = stoichiometry.frac_mole_to_weight(np.asarray(list(nfrac.values())), MM)
        return dict(zip(nfrac.keys(), wfrac))

    def elemental_molefractions(self):
        return self.molefractions()

    def elemental_equivalents(self):
        return self.equivalents()

    def elemental_massfractions(self):
        return self.massfractions()

    def fractions(self, fractype):
        if fractype == types.fraction.mole:
            return self.molefractions()
        elif fractype == types.fraction.volume:
            raise ValueError("Cannot calculate elemental volume fractions")
        else:
            return self.massfractions()

    @property
    def nelements(self):
        return len(self._elements)

    def markabsorber(self, symb=None, shells=None, fluolines=None, energybounds=None):
        """
        Args:
            symb(str): element symbol
        """
        for e in self._elements:
            if energybounds is None:
                energybounds2 = None
            else:
                energybounds2 = energybounds
            e.markabsorber(
                symb=symb,
                shells=shells,
                fluolines=fluolines,
                energybounds=energybounds2,
            )

    def unmarkabsorber(self):
        for e in self._elements:
            e.unmarkabsorber()

    def hasabsorbers(self):
        return any([e.isabsorber for e in self._elements])

    def markscatterer(self, name=None):
        """
        Args:
            name(str): compound name
        """
        if name is None:
            self.isscatterer = True
        else:
            self.isscatterer = self == name

    def unmarkscatterer(self):
        self.isscatterer = False

    def markinfo(self):
        yield "{}".format(self.name)
        yield " Scatterer: {}".format("yes" if self.isscatterer else "no")
        for e in self._elements:
            for s in e.markinfo():
                yield " {}".format(s)

    def _crosssection(self, method, E, fine=False, decomposed=False, **kwargs):
        """Calculate compound cross-sections"""
        if self._cs_scattering(method) and not self.isscatterer:
            if decomposed:
                return {}
            else:
                return np.zeros_like(E, dtype=float)
        e_wfrac = self.massfractions()
        if hasattr(self, "structure") and fine:
            environ = self
        else:
            environ = None
        kwargs["environ"] = environ
        if decomposed:
            ret = {}
            for e, w in e_wfrac.items():
                cs = getattr(e, method)(E, **kwargs)
                ret[e] = {"w": w, "cs": cs}
        else:
            if self._cs_dict(method):
                ret = {}
                for e, w in e_wfrac.items():
                    cs = getattr(e, method)(E, **kwargs)
                    if not cs:
                        continue
                    for k, v in cs.items():
                        ret[k] = ret.get(k, 0) + w * v
            else:
                ret = sum(
                    w * getattr(e, method)(E, **kwargs) for e, w in e_wfrac.items()
                )
        return ret

    def get_energy(self, energyrange, defaultinc=1):
        """Get absolute energies (keV) from a relative energy range (eV)

        Args:
            energyrange(np.array): energy range relative to the absorber's edges (if any)
        """
        for e in self._elements:
            ret = e.get_energy(energyrange, defaultinc=defaultinc)
            if ret is not None:
                return ret
        return None

    def topymca(self, cfg, defaultthickness=1e-4):
        r = self.massfractions()
        massfractions = list(r.values())
        names = ["{}1".format(e) for e in r]
        matname = self.pymcaname
        cfg["materials"][matname] = {
            "Comment": self.pymcacomment,
            "CompoundFraction": massfractions,
            "Thickness": defaultthickness,
            "Density": self.density,
            "CompoundList": names,
        }
        return matname

    def tofisx(self, cfg, defaultthickness=1e-4):
        r = self.massfractions()
        massfractions = list(r.values())
        names = ["{}1".format(e) for e in r]
        matname = self.pymcaname
        o = fisx.Material(matname, self.density, defaultthickness, self.pymcacomment)
        o.setCompositionFromLists(names, massfractions)
        cfg.addMaterial(o, errorOnReplace=False)
        return matname

    def tocompound(self, name=None):
        return self

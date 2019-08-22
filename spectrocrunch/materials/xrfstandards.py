# -*- coding: utf-8 -*-

from . import compoundfromname
from . import compoundfromlist
from . import mixture
from . import element
from . import multilayer
from . import types
from . import pymca
from ..utils.classfactory import with_metaclass
from ..utils import instance

import numpy as np


class Standard(with_metaclass(multilayer.Multilayer)):

    def __init__(self, extraelements=None, **kwargs):
        if extraelements is None:
            self.extraelements = None
        else:
            if instance.isstring(extraelements):
                extraelements = [extraelements]
            self.extraelements = [element.Element(el) for el in extraelements]
        super(Standard, self).__init__(**kwargs)

    def addtopymca(self, setup, cfg):
        super(Standard, self).addtopymca(setup, cfg)
        if self.extraelements is not None:
            self.addtopymca_shells(setup, cfg, self.extraelements)


class ThinFilmStandard(Standard):

    def __init__(self, arealdensity, substrate, substratethickness,
                 name=None, filmthickness=None, **kwargs):
        arealdensity_film = dict(
            zip(map(element.Element, arealdensity.keys()), arealdensity.values()))

        if filmthickness is None:
            vacuum = compoundfromname.compoundfromname('dummy')
            wfrac_substrate = substrate.elemental_massfractions()

            elfilm = set(arealdensity_film.keys())
            elsubstrate = set(wfrac_substrate.keys())
            if elfilm & elsubstrate:
                raise RuntimeError(
                    "Film and substrate cannot contain the same elements")

            #kwargs["material"] = substrate
            kwargs["material"] = vacuum
            kwargs["thickness"] = substratethickness
            if "extraelements" not in kwargs:
                kwargs["extraelements"] = []
            kwargs["extraelements"].extend(elfilm)
            kwargs["extraelements"].extend(elsubstrate)

            self._arealdensities = arealdensity_film
            m = substrate.density*substratethickness
            self._arealdensities.update(
                dict(zip(wfrac_substrate.keys(), np.array(wfrac_substrate.values())*m)))
        else:
            ad = np.array(arealdensity_film.values())
            sad = ad.sum()
            density = sad/filmthickness
            massfractions = ad/sad
            film = compoundfromlist.CompoundFromList(arealdensity_film.keys(), massfractions,
                                                     types.fraction.mass, density=density, name=name)
            kwargs["material"] = [film, substrate]
            kwargs["thickness"] = [filmthickness, substratethickness]

            self._arealdensities = None

        super(ThinFilmStandard, self).__init__(**kwargs)

    def arealdensity(self):
        if self._arealdensities is None:
            return super(ThinFilmStandard, self).arealdensity()
        else:
            return self._arealdensities


class InfThickStandard(Standard):
    pass


class AXOID21_1(ThinFilmStandard):
    aliases = ["RF7-200-S2371-03", "id21_room"]

    def __init__(self, **kwargs):
        name = "RF7-200-S2371-03"

        # ng/mm^2 -> g/cm^2
        arealdensity = {"Pb": 7.7e-7,
                        "La": 9.0e-7,
                        "Pd": 1.9e-7,
                        "Mo": 0.9e-7,
                        "Cu": 2.4e-7,
                        "Fe": 4.0e-7,
                        "Ca": 11.4e-7}
        filmthickness = kwargs.pop("filmthickness", None)  # cm
        substrate = compoundfromname.compoundfromname("silicon nitride")
        substratethickness = 200e-7  # cm

        ultralene = compoundfromname.compoundfromname("ultralene")
        attenuators = [["SampleCover", ultralene, 4.064e-4],
                       ["BeamFilter0", ultralene, 4.064e-4]]
        for k in attenuators:
            kwargs["geometry"].addattenuator(*k)

        super(AXOID21_1, self).__init__(arealdensity, substrate,
                                        substratethickness, name=name, filmthickness=filmthickness, **kwargs)


class AXOID21_2(ThinFilmStandard):
    aliases = ["RF8-200-S2454-17", "id21_cryo"]

    def __init__(self, **kwargs):
        name = "RF8-200-S2454-17"

        # ng/mm^2 -> g/cm^2
        arealdensity = {"Pb": 6.3e-7,
                        "La": 7.6e-7,
                        "Pd": 2.3e-7,
                        "Mo": 0.7e-7,
                        "Cu": 2.6e-7,
                        "Fe": 4.1e-7,
                        "Ca": 25.1e-7}
        filmthickness = kwargs.pop("filmthickness", None)  # cm
        substrate = compoundfromname.compoundfromname("silicon nitride")
        substratethickness = 200e-7  # cm

        ultralene = compoundfromname.compoundfromname("ultralene")
        attenuators = [["SampleCover", ultralene, 4.064e-4],
                       ["BeamFilter0", ultralene, 4.064e-4]]
        for k in attenuators:
            kwargs["geometry"].addattenuator(*k)

        super(AXOID21_2, self).__init__(arealdensity, substrate,
                                        substratethickness, name=name, filmthickness=filmthickness, **kwargs)


class AXOID16b_1(ThinFilmStandard):
    aliases = ["RF8-200-S2453"]

    def __init__(self, **kwargs):
        name = "RF8-200-S2453"

        # ng/mm^2 -> g/cm^2
        arealdensity = {"Pb": 5.9e-7,
                        "La": 10.3e-7,
                        "Pd": 1.2e-7,
                        "Mo": 0.7e-7,
                        "Cu": 2.4e-7,
                        "Fe": 3.9e-7,
                        "Ca": 20.3e-7}
        filmthickness = kwargs.pop("filmthickness", None)  # cm
        substrate = compoundfromname.compoundfromname("silicon nitride")
        substratethickness = 200e-7  # cm

        super(AXOID16b_1, self).__init__(arealdensity, substrate,
                                         substratethickness, name=name, filmthickness=filmthickness, **kwargs)


class NIST612(Standard):

    def __init__(self, **kwargs):
        major = {compoundfromname.compoundfromname("silica"): 0.72,
                 compoundfromname.compoundfromname("sodium oxide"): 0.14,
                 compoundfromname.compoundfromname("lime"): 0.12,
                 compoundfromname.compoundfromname("corundum"): 0.02}

        glass = mixture.Mixture(major.keys(), major.values(
        ), types.fraction.mass, name="glass").tocompound()

        certified = {"Sb": 34.9, "As": 37.4, "Ba": 38.6, "Cd": 29.9, "Cr": 35.0, "Fe": 51., "Pb": 38.57, "Mn": 37.7,
                     "Ni": 38.8, "Rb": 31.4, "Se": 16.1, "Ag": 22.0, "Sr": 78.4, "Th": 37.79, "U": 37.38}

        reference = {"Co": 35.5, "Cu": 37.7, "Ta": 15.7, "Ti": 50.1}

        info = {"B": 32., "Ce": 39., "Dy": 35., "Er": 39., "Eu": 36., "Gd": 39., "Au": 5., "La": 36., "Li": 40., "Nd": 36.,
                "K": 64., "Sm": 39., "Yb": 42.}

        elements = certified.keys()+reference.keys()+info.keys()
        fractions = certified.values()+reference.values()+info.values()
        glass.addelements(elements, np.array(
            fractions)/1e6, types.fraction.mass)

        if "extraelements" not in kwargs:
            kwargs["extraelements"] = []
        # kwargs["extraelements"].extend(reference.keys()+info.keys())

        super(NIST612, self).__init__(
            material=glass, thickness=300e-4, **kwargs)


factory = Standard.factory
registry = Standard.clsregistry

from ...geometries import xrf as xrfgeometries
from ...sources import xray as xraysources
from ...detectors import xrf as xrfdetectors
from .. import compoundfromname
from .. import compoundfromformula
from .. import element
from .. import mixture
from .. import types
from .. import multilayer
from .. import pymca


def geometry():
    source = xraysources.factory("synchrotron")
    detector = xrfdetectors.factory("leia")
    return xrfgeometries.factory(
        "sxm120", detector=detector, source=source, detectorposition=-15.0
    )


def single(**kwargs):
    # Cover elements, compounds and mixtures
    # Not too many lines for speed
    ca = element.Element("Ca")
    calcite = compoundfromname.compoundfromname("calcite")
    hematite = compoundfromname.compoundfromname("hematite")
    goethite = compoundfromname.compoundfromname("goethite")
    mix = mixture.Mixture(
        [ca, hematite, goethite, calcite],
        [0.25, 0.25, 0.25, 0.25],
        types.fraction.mass,
        name="iron oxides",
    ).tocompound()  # TODO: fix pymca bug
    sample = multilayer.Multilayer(mix, 1e-4, geometry=geometry())
    kwargs["energy"] = kwargs.get("energy", 8.0)
    pymcahandle = pymca.PymcaHandle(sample=sample, **kwargs)
    return pymcahandle


def simple(**kwargs):
    # Cover elements, compounds and mixtures
    # Not too many lines for speed
    ca = element.Element("Ca")
    calcite = compoundfromname.compoundfromname("calcite")
    hematite = compoundfromname.compoundfromname("hematite")
    goethite = compoundfromname.compoundfromname("goethite")
    mix = mixture.Mixture(
        [hematite, goethite], [0.5, 0.5], types.fraction.mass, name="iron oxides"
    ).tocompound()  # TODO: fix pymca bug
    sample = multilayer.Multilayer(
        [ca, mix, calcite], [2e-5, 7e-5, 10e-5], geometry=geometry()
    )
    kwargs["energy"] = kwargs.get("energy", 8.0)
    pymcahandle = pymca.PymcaHandle(sample=sample, **kwargs)
    return pymcahandle


def complex(**kwargs):
    # Cover L and M lines
    ca = element.Element("Ca")
    calcite = compoundfromname.compoundfromname("calcite")
    hematite = compoundfromname.compoundfromname("hematite")
    goethite = compoundfromname.compoundfromname("goethite")
    mix = mixture.Mixture(
        [hematite, goethite], [0.5, 0.5], types.fraction.mass, name="iron oxides"
    ).tocompound()  # TODO: fix pymca bug
    compound1 = compoundfromformula.CompoundFromFormula("PbCe", density=6.0)
    sample = multilayer.Multilayer(
        [ca, mix, compound1, calcite], [2e-5, 3e-5, 1e-5, 10e-5], geometry=geometry()
    )
    kwargs["energy"] = kwargs.get("energy", 8.0)
    pymcahandle = pymca.PymcaHandle(sample=sample, **kwargs)
    return pymcahandle

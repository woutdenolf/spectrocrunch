import numpy as np
from . import xrmc
from ..sources import xray as xraysources
from ..materials import pymca
from ..materials import multilayer
from ..detectors import xrf as xrfdetectors
from ..geometries import xrf as xrfgeometries
from ..materials import compoundfromformula
from ..materials import compoundfromdb
from ..materials import compoundfromlist
from ..materials import mixture
from ..materials import types
from ..materials import xrfstandards
from ..io import localfs


def define_source(source=None, energy=7.5, flux=1e10, time=0.1, **kwargs):
    if source is None:
        source = xraysources.factory("synchrotron")
    pymcahandle = pymca.PymcaHandle(
        energy=energy,
        flux=flux,
        time=time,
        linear=True,
        escape=True,
        snip=True,
        scatter=True,
        **kwargs,
    )
    return source, pymcahandle


def world_addsource(world, pymcahandle, beamsize=None, fullfield=False):
    distance = 42e2  # doesn't really matter
    # Beamsize at the sample
    if beamsize is None:
        if fullfield:
            beamsize = 100e-4  # cm
        else:
            beamsize = 1e-4  # cm
    world.define_source(
        flux=pymcahandle.flux.to("Hz").magnitude,
        energy=pymcahandle.energy,
        distance=distance,
        beamsize=beamsize,
    )


def define_xrfgeometry(source):
    detector = xrfdetectors.factory("leia")
    geometry = xrfgeometries.factory("sxm120", detector=detector, source=source)
    return geometry


def world_addsdd(
    world,
    pymcahandle,
    distance=None,
    polar=None,
    azimuth=None,
    activearea=None,
    hoffset=0,
    voffset=0,
    convoluted=True,
    pulseproctime=0,
    multiplicity=10,
    emin=None,
    emax=None,
):
    forcedetect = True
    poissonnoise = False
    geometry = pymcahandle.sample.geometry
    ebinsize = geometry.detector.mcagain  # keV
    if distance is None:
        distance = geometry.distance.to("cm").magnitude
    if activearea is None:
        activearea = geometry.detector.activearea.to("cm**2").magnitude
    if polar is None:
        polar = geometry.scatteringangle
    if azimuth is None:
        azimuth = geometry.detector_azimuth
    if convoluted:
        response = {
            "material": geometry.detector.material,
            "thickness": geometry.detector.thickness,
            "noise": geometry.detector.mcanoise,
            "fano": geometry.detector.mcafano,
            "pulseproctime": pulseproctime,
        }
    else:
        response = {}
    world.add_xrfdetector(
        distance=distance,
        activearea=activearea,
        polar=polar,
        azimuth=azimuth,
        hoffset=hoffset,
        voffset=voffset,
        emin=emin,
        emax=emax,
        ebinsize=ebinsize,
        poissonnoise=poissonnoise,
        forcedetect=forcedetect,
        multiplicity=multiplicity,
        response=response,
        time=pymcahandle.time.to("s").magnitude,
    )
    if convoluted:
        geometry.detector.mcagain = world.detector.mcagain
        geometry.detector.mcazero = world.detector.mcazero
    dx, dy = world.detector.pixelsize
    ox = 0
    oy = 0
    for layer in geometry.detectorfilters(include_atmosphere=False):
        world.detector.add_layer(
            material=layer["material"],
            thickness=-layer["thickness"],
            dhor=dx,
            dvert=dy,
            ohor=ox,
            overt=oy,
        )


def define_sample(geometry, name=None):
    if name == "cell":
        cell = compoundfromdb.compoundfromdb("ecoli dry")
        substrate = compoundfromdb.compoundfromdb("silicon nitride")
        # Dope cell with elements (1e-6 = 1 ppm)
        cell.addelements(["Ce"], [100e-6], types.fraction.mole)
        # Add sample cover
        ultralene = compoundfromdb.compoundfromdb("ultralene")
        attenuators = [
            ["SampleCover", ultralene, 4.064e-4],
            ["BeamFilter0", ultralene, 4.064e-4],
        ]
        for k in attenuators:
            geometry.addattenuator(*k)
        # Multilayer of cell on substrate
        sample = multilayer.Multilayer(
            material=[cell, substrate],
            thickness=[10e-4, 0.2e-4],
            geometry=geometry,
            name="cell",
        )
    elif name == "plant":
        plant = compoundfromdb.compoundfromdb("Cellulose Nitrate")
        # Dope cell with elements (1e-6 = 1 ppm)
        elements = ["Ce", "P", "S", "Cl", "K", "Ca", "Mn", "Fe"]
        mfracs = [200e-6, 2000e-6, 1000e-6, 100e-6, 10000e-6, 5000e-6, 50e-6, 100e-6]
        plant.addelements(elements, mfracs, types.fraction.mole)
        # Add sample cover
        ultralene = compoundfromdb.compoundfromdb("ultralene")
        attenuators = [
            ["SampleCover", ultralene, 4.064e-4],
            ["BeamFilter0", ultralene, 4.064e-4],
        ]
        for k in attenuators:
            geometry.addattenuator(*k)
        # Multilayer of cell on substrate
        sample = multilayer.Multilayer(
            material=[plant, ultralene],
            thickness=[20e-4, 4.064e-4],
            geometry=geometry,
            name="plant",
        )
    elif name == "imaging":
        c1 = compoundfromdb.compoundfromdb("hematite")
        c2 = compoundfromformula.CompoundFromFormula("PbSO4", density=6.29)
        c3 = compoundfromlist.CompoundFromList(
            ["Ca", "C", "O"],
            [1, 1, 3],
            types.fraction.mole,
            density=2.71,
            name="calcite",
        )
        toplayer = c1
        substrate = mixture.Mixture(
            [c2, c3], [1, 1], types.fraction.mass, name="Substrate"
        )
        sample = multilayer.Multilayer(
            material=[toplayer, substrate],
            thickness=[10e-4, 20e-4],
            geometry=geometry,
            name="imaging",
        )
    elif name == "axo":
        sample = xrfstandards.factory(
            "id21_room",
            geometry=geometry,
            filmthickness=1e-7,
            extraelements=["Ni", "Cl", "Cr", "Al"],
        )
    else:
        c1 = compoundfromdb.compoundfromdb("hematite")
        c2 = compoundfromformula.CompoundFromFormula("PbSO4", density=6.29)
        c3 = compoundfromlist.CompoundFromList(
            ["Ca", "C", "O"],
            [1, 1, 3],
            types.fraction.mole,
            density=2.71,
            name="calcite",
        )
        homogen = mixture.Mixture(
            [c1, c2, c3], [1, 1, 1], types.fraction.mass, name="mixture"
        )
        sample = multilayer.Multilayer(
            material=[homogen], thickness=[30e-4], geometry=geometry, name="simple"
        )
    return sample


def world_addsample(world, sample, fullfield=False):
    projbeamsizesample = 2 * (world.source.distance) * np.tan(world.source.divergence)
    assert world.source.beamsize == projbeamsizesample
    if fullfield:
        # Make sure the entire sample fits in the beam footprint
        samplesize = projbeamsizesample / 2.0
    else:
        # Make sure the sample is larger than the beam footprint
        samplesize = max(projbeamsizesample * 10, 0.5)

    # All layers the same size and beam throught sample center
    nlayers = len(sample)
    dxs = [samplesize] * nlayers
    dys = [samplesize] * nlayers
    oxs = [0] * nlayers
    oys = [0] * nlayers

    if sample.name == "imaging":
        # Make first layer smaller
        dxs[0] = 3 * samplesize / 4.0
        dys[0] = samplesize / 4.0
        oxs[0] = (samplesize - dxs[0]) / 2.0
        oys[0] = 0  # -(samplesize-dys[0])/2.

    for layer, dx, dy, ox, oy in zip(sample, dxs, dys, oxs, oys):
        world.addlayer(
            material=layer.material,
            thickness=layer.thickness,
            dhor=dx,
            dvert=dy,
            ohor=ox,
            overt=oy,
        )
    return sample


def run(world, interactions, simulate=True, plot=True, ylog=False):
    path = localfs.Path(world.main.path)
    if simulate:
        path.remove(recursive=True)
        world.finalize(interactions=interactions)
        if not world.simulate():
            # path.ls(recursive=True)
            raise RuntimeError("Simulation failed")
    else:
        world.finalize(interactions=interactions)
    # if plot:
    #    path.ls(recursive=True)
    if world.detector.TYPE == "detectorconvolute":
        basename = world.detector.name + "_convoluted"
    else:
        basename = world.detector.name + "_lines"
    data, info = xrmc.loadxrmcresult(world.detector.outpath, basename)
    if plot:
        # info.pop('xenergy') # in case you want channels
        xrmc.showxrmcresult(data, ylog=ylog, **info)
    return data, info

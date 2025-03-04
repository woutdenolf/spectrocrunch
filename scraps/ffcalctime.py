import os

from spectrocrunch.simulation import calcnoise
from spectrocrunch.materials import multilayer
from spectrocrunch.materials.compoundfromname import compoundfromname
from spectrocrunch.materials.mixture import Mixture

from spectrocrunch.fullfield import create_hdf5_imagestacks as ff

import numpy as np
from uncertainties import unumpy
import matplotlib.pyplot as plt


def getdata(sourcepath, radix, roi=((0, None), (0, None))):
    config = {
        # EDF header
        "frametimelabel": "exposure_time",
        "frametimedefault": 0.0,
        "roilabel": "roi",
        "nbframeslabel": "nb_frames",
        "stacklabel": "energy",
        # Defaults
        "dtype": "np.float32",  # must be floating point!
        "darkcurrentzero": 90.0,
        "darkcurrentgain": 1.0,
        "stackdim": 2,
        # Data
        "darklist": map(
            lambda xy: os.path.join(xy[0], xy[1] + "*_dark_*.edf"),
            zip(sourcepath, radix),
        ),
        "datalist": map(
            lambda xy: os.path.join(xy[0], xy[1] + "*_data_*.edf"),
            zip(sourcepath, radix),
        ),
        "flatlist": map(
            lambda xy: os.path.join(xy[0], xy[1] + "*_ref_*.edf"),
            zip(sourcepath, radix),
        ),
        "beforeafter": True,
        # Output
        "roi": roi,
        "rebin": [1, 1],
    }

    # Dark, data and flat libraries
    darklib = ff.darklibrary(config)
    data, flat1, flat2, keyindices, stackaxes = ff.dataflatlibrary(config)
    keys = data.keys()

    N = np.empty(len(stackaxes[2]["data"]), dtype=np.float32)
    D = np.empty(len(stackaxes[2]["data"]), dtype=np.float32)
    N0 = np.empty(len(stackaxes[2]["data"]), dtype=np.float32)
    D0 = np.empty(len(stackaxes[2]["data"]), dtype=np.float32)
    energy = stackaxes[2]["data"]

    for i in range(len(keyindices)):
        key = keys[keyindices[i]]
        img, tframe_data, nframes_data, dark_data = ff.getsingleimage(
            data[key][0], darklib, config
        )
        nframes_dark = dark_data["nframes"]
        dark_data = dark_data["data"] / dark_data["nframes"] * nframes_data

        # Normalize
        flat, tframe_flat, nframes_flat, dark_flat = ff.getsingleimage(
            flat1[key][0], darklib, config
        )
        if flat2 is not None:
            tmp, _, nframes_flat2, _ = ff.getsingleimage(flat2[key][0], darklib, config)
            flat += tmp
            nframes_flat += nframes_flat2
        dark_flat = dark_flat["data"] / dark_flat["nframes"] * nframes_flat

        N[i] = img.mean()
        N0[i] = flat.mean()
        D[i] = dark_data.mean()
        D0[i] = dark_flat.mean()

    return (
        energy,
        N,
        N0,
        D,
        D0,
        tframe_data,
        nframes_data,
        tframe_flat,
        nframes_flat,
        nframes_dark,
    )


def hg107():

    ultralene = compoundfromname("ultralene")
    tape = compoundfromname("kapton")
    prussianblue = compoundfromname("prussianblue")
    hydrocerussite = compoundfromname("hydrocerussite")

    _ = multilayer.Multilayer(
        material=[ultralene, prussianblue, hydrocerussite, tape],
        thickness=[4.0, 10.0, 10.0, 5.0],
    )

    # getsignal(sourcepath, radix, indices)


def hg64():
    # Estimate sample
    ultralene = compoundfromname("ultralene")
    tape = compoundfromname("kapton")
    cadmiumsulfide = compoundfromname("cadmium sulfide")
    sample = multilayer.Multilayer(
        material=[ultralene, cadmiumsulfide, tape], thickness=[4.0, 10.0, 5.0]
    )

    # Experimental data
    sourcepath = [
        "/data/id21/store/backup_visitor/2016/hg64/id21/7913_50RH_ffCd/ff/map1"
    ]
    radix = ["map1"]
    x = 966
    y = 310
    d = 5
    (
        energy,
        N,
        N0,
        D,
        D0,
        tframe_data,
        nframes_data,
        tframe_flat,
        nframes_flat,
        nframes_dark,
    ) = getdata(sourcepath, radix, roi=((y - d, y + d + 1), (x - d, x + d + 1)))
    xanes = -np.log((N - D) / (N0 - D0))

    plt.figure()
    plt.plot(energy, xanes)
    plt.xlabel("Energy (keV)")
    plt.ylabel("Absorbance")

    # Subset of experimental data
    energysel = np.array([3.49, 3.5, 3.68, 3.72])
    xanessel = energysel * 0
    for i in range(len(energysel)):
        j = np.argmin(np.abs(energy - energysel[i]))
        energysel[i] = energy[j]
        xanessel[i] = xanes[j]

    # For layer thicknesses on subset
    sample.refinethickness(energysel, xanessel, layerfixed=[0])
    for layer in range(sample.nlayers):
        print("{}: {} um".format(sample.material[layer], sample.thickness[layer]))

    plt.plot(energysel, xanessel, "o")
    energy = np.linspace(3.48, 3.8, 100)
    xanes = sample.absorbance(energy)
    plt.plot(energy, xanes)

    # Noise propagation
    I0 = 1e0

    print(tframe_data, nframes_data, tframe_flat, nframes_flat, nframes_dark)
    N, N0, D, D0 = calcnoise.id21_ffnoise(
        I0,
        energy,
        sample,
        tframe_data,
        nframes_data,
        tframe_flat,
        nframes_flat,
        nframes_dark,
    )

    XAS = -unumpy.log((N - D) / (N0 - D0))

    signal = unumpy.nominal_values(XAS)
    # noise = unumpy.std_devs(XAS)
    plt.plot(energy, signal)

    plt.show()


def hg94():
    ultralene = compoundfromname("ultralene")
    binder = compoundfromname("linseed oil")
    pigment = compoundfromname("lazurite")
    mix = Mixture([binder, pigment], [0.5, 0.5])

    sample = multilayer.Multilayer(
        material=[ultralene, mix, ultralene],
        thickness=[4.064e-4, 10e-4, 4.064e-4],
        fixed=[True, False, True],
    )

    flux = 1e5
    energy = np.linspace(2.46, 2.65, 140)
    params = {
        "tframe_data": 0.07,
        "tframe_flat": 0.07,
        "nframe_data": 100,
        "nframe_flat": 100,
        "nframe_dark": 10,
    }

    ffsetup = calcnoise.id21_ffsetup(sample=sample)

    signal, noise = ffsetup.xanes(flux, energy, **params)
    signal = np.random.normal(signal, noise)
    _ = plt.plot(energy, signal)
    plt.xlabel("Energy (keV)")
    plt.ylabel("Absorbance")
    plt.show()


if __name__ == "__main__":

    hg94()

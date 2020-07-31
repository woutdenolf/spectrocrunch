# -*- coding: utf-8 -*-

import os, sys

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spectrocrunch.io.spec import spec
from spectrocrunch.materials.compoundfromformula import compoundfromformula as compound

import matplotlib.pyplot as plt
import numpy as np
import xraylib
import scipy.optimize


def errorf_attzscan(p, att, iodet, idet, time):
    I0, b0, g0, b, g = p
    ret = (b0 + g0 * I0 * time * att - iodet) / iodet
    ret = np.append(ret, (b + g * I0 * att - idet) / idet)
    return ret


def attzscan(specfile, scannumber, outpath, withsample=False, iodettype=2):
    """
        I0: photons per second
        Photons: # photons detected by idet
        attfoil = exp(-mu.rho.d)
        attiodet = exp(-mu.rho.d)
        attsample = exp(-mu.rho.d)

        Photons = I0*time*attfoil*attiodet*attsample

        idet  = b  +  g*Photons/time
              = b  +  g*I0*attfoil*attiodet*attsample

        I0 = Photons/(time*attfoil*attiodet*attsample)

        iodet = b0 + g0*I0*attfoil
              = b0 + g0*Photons/(time*attiodet*attsample)
        
    """

    # Data
    f = spec(specfile)
    result = f.getdata2(
        scannumber, ["attz", "iodet", "idet", "Seconds", "Photons", "Storage Ring A"]
    )
    energy = np.array(f.getmotorvalues(scannumber, ["Energy MONO"]))
    ind = np.argsort(result[:, 0])
    attz = result[ind, 0]
    iodet = result[ind, 1]
    idet = result[ind, 2]
    time = result[ind, 3]
    photons = result[ind, 4]
    srcur = result[ind, 5]

    # Interpolated data
    n = 11
    attzi = np.arange(n) * 6 + 1.5
    attzi[0] -= 1
    attzi = attzi[attzi > attz[0]]
    attzi = attzi[attzi < attz[-1]]
    attzi = attzi[attzi < 55]
    n = len(attzi)
    iodeti = np.interp(attzi, attz, iodet)
    ideti = np.interp(attzi, attz, idet)
    timei = np.interp(attzi, attz, time)
    photonsi = np.interp(attzi, attz, photons)
    srcuri = np.interp(attzi, attz, srcur)

    # Attenuation
    Al = compound("Al", 0, name="Al")
    Ti = compound("Ti", 0, name="Ti")
    Si3N4 = compound("Si3N4", 3.44, name="Si3N4")
    ultralene = compound("C22H10N205", 1.42, "Ultralene")

    thickness = np.array([0, 6, 12, 20, 30, 40, 50, 100, 150, 25, 50]) * 1e-4
    comp = [Al, Al, Al, Al, Al, Al, Al, Al, Al, Ti, Ti]
    attfoil = np.empty(n)
    for i in range(n):
        attfoil[i] = np.exp(
            -comp[i].mass_att_coeff(energy)[0] * comp[i].density * thickness[i]
        )
    attiodet = np.exp(-Si3N4.mass_att_coeff(energy)[0] * Si3N4.density * (0.5 * 1e-4))
    if iodettype == 1:
        attiodet = np.exp(-Ti.mass_att_coeff(energy)[0] * Ti.density * (1e-4))
    if withsample:
        attsample = np.exp(
            -ultralene.mass_att_coeff(energy)[0] * ultralene.density * (4 * 1e-4)
        )
    else:
        attsample = 1

    # Gain and offset of idet
    A = np.vstack([photons, np.ones(len(photons))]).T
    g, b = np.linalg.lstsq(A, idet)[0]
    A = np.vstack([photons / (attiodet * attsample), np.ones(len(photons))]).T
    g0, b0 = np.linalg.lstsq(A, iodet)[0]
    # b = 0
    # b0 = 200

    # Flux based on the flux w.o. sample and the SR current
    I0th = photonsi / (attfoil * attiodet * attsample)
    I0 = I0th[0] / (srcuri[0] * timei[0]) * srcuri * timei[0]

    # Plot
    print("Scan: {}".format(scannumber))
    print("Cur: {} mA".format(srcuri[0]))
    print("Energy: {} keV".format(energy[0]))
    print("I0: {:e} ph".format(I0[0]))
    print("b0: {:e} cts".format(b0))
    print("g0: {:e} cts/ph".format(g0))
    print("b: {:e} cts".format(b))
    print("g: {:e} cts/ph".format(g))
    print("")

    fig = plt.figure(1)
    plt.title("{}".format(specfile))
    plt.plot(I0th / srcuri / timei, label="{}".format(scannumber))
    plt.xlabel("attz")
    plt.ylabel("I0/SR/$\Delta$t (ph/s/mA)")
    plt.legend(loc=3)
    fig.savefig(os.path.join(outpath, "I0.png"), bbox_inches="tight", dpi=300)

    fig = plt.figure(2)
    plt.title("{}".format(specfile))
    p = plt.plot(photons, idet, "o", label="{}".format(scannumber))
    plt.plot(photons, b + g * photons, color=p[0].get_color())
    plt.xlabel("Photons")
    plt.ylabel("idet")
    plt.legend(loc=2)
    fig.savefig(os.path.join(outpath, "idet.png"), bbox_inches="tight", dpi=300)

    fig = plt.figure(3)
    plt.title("{}".format(specfile))
    p = plt.plot(photons, iodet, "o", label="{}".format(scannumber))
    plt.plot(
        photons, b0 + g0 * photons / (attiodet * attsample), color=p[0].get_color()
    )
    plt.xlabel("Photons")
    plt.ylabel("iodet")
    plt.legend(loc=2)
    fig.savefig(os.path.join(outpath, "iodet.png"), bbox_inches="tight", dpi=300)

    fig = plt.figure()
    p = plt.plot(attz, iodet, label="iodet")
    plt.plot(
        attzi, b0 + g0 * I0 * attfoil, "o", label="iodet(I0)", color=p[0].get_color()
    )
    p = plt.plot(attz, idet, label="idet")
    plt.plot(
        attzi,
        b + g * I0 * attfoil * attiodet * attsample,
        "o",
        label="idet(I0)",
        color=p[0].get_color(),
    )
    plt.xlabel("attz")
    plt.title("{}: {}".format(specfile, scannumber))
    plt.legend()
    fig.savefig(
        os.path.join(outpath, "{}.png".format(scannumber)), bbox_inches="tight", dpi=300
    )

    return [scannumber, b0, g0, b, g]


def idet(v):
    """
        ph_E: energy in keV
        ph_I: idet * pow(10, -5-ph_gain)
        ph_coeff: linear att. coeff. air
        
        photons/s = ph_I * 1e3 / (ph_E * 1.6E-16 * exp(-ph_coeff * PH_DIST) * ph_factor * ph_PTB  )
    """
    pass


if __name__ == "__main__":
    specfile = "/data/id21/inhouse/16mar/BLC_quanti/data/16033001.dat"
    outpath = "/data/id21/inhouse/16mar/BLC_quanti/results/"

    c = []
    c.append(attzscan(specfile, 111, outpath, withsample=True))
    c.append(attzscan(specfile, 112, outpath, withsample=False))
    c.append(attzscan(specfile, 121, outpath, withsample=False))
    c.append(attzscan(specfile, 122, outpath, withsample=True))
    c.append(attzscan(specfile, 140, outpath, withsample=False))
    c.append(attzscan(specfile, 141, outpath, withsample=False))
    c = np.array(c)
    np.savetxt(os.path.join(outpath, "gains.csv"), c, delimiter=";", fmt="%s")

    # plt.show()

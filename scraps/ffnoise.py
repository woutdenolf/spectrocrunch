# -*- coding: utf-8 -*-

import os, sys

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spectrocrunch.io.spec import spec
import matplotlib.pyplot as plt
import numpy as np

from spectrocrunch.materials.compoundfromformula import compoundfromformula as compound

# from spectrocrunch.materials.compoundfromcif import compoundfromcif as compound
from spectrocrunch.materials.mixture import mixture
from spectrocrunch.materials.types import fraction

import xraylib
import silx.math.fit as fit


def gettransmissionxanes(filename, scannumber):
    f = spec(filename)
    result = f.getdata2(scannumber, ["arr_energyM", "arr_iodet", "arr_idet"])
    energy = result[:, 0]
    iodet = result[:, 1]
    idet = result[:, 2]
    return energy, iodet, idet


def muL_combination(x, thickness, w1):
    energy, b, rho = x
    w2 = 1 - w1
    rho = 1 / (w1 / rho[0] + w2 / rho[1])
    return thickness * rho * (w1 * b[0] + w2 * b[1])


def getcalculated(energy, e1, e2, e3, e4, absorbance_measured=None):

    compound1 = compound("Na4Ca4Al6Si6O32S2", 2.4, name="pigment")
    compound2 = compound("C22H10N205", 1.43, name="kapton")
    compounds = [compound1, compound2]
    rho = [c.density for c in compounds]
    mask = (energy >= e1) & (energy <= e2)
    mask |= (energy >= e3) & (energy <= e4)
    b = [c.mass_att_coeff(energy[mask]) for c in compounds]
    print(energy[mask])

    thickness = 5e-4
    w1 = 1.0

    if absorbance_measured is not None:
        x = energy[mask], b, rho
        y = absorbance_measured[mask]
        ind = y < 0
        y[ind] = 0
        ysigma = np.sqrt(y)
        ysigma[ind] = 1

        p0 = (thickness, w1)
        # constraints = [[fit.CFREE,0,0],[fit.CQUOTED,0,1]]
        constraints = [[fit.CFREE, 0, 0], [fit.CFIXED, 0, 0]]

        p, cov_matrix, info = fit.leastsq(
            muL_combination,
            x,
            y,
            p0,
            sigma=ysigma,
            constraints=constraints,
            full_output=True,
        )
        var = np.diag(cov_matrix)

        S = np.diag(1 / np.sqrt(var))
        cor_matrix = S.dot(cov_matrix).dot(S)
        print("R(d,w1) = {}".format(cor_matrix[0, 1]))

        thickness, w1 = p
        thicknessvar, w1var = var
        w2 = 1 - w1
        w2var = w1var

    m = mixture([compound1, compound2], [w1, w2], fraction.mass)
    muL = m.mass_att_coeff(energy, decomposed=False) * m.density() * thickness

    print("Density = {} g/cm^3".format(m.density()))
    print(
        "Thickness = {} +/- {} um".format(thickness * 1e4, np.sqrt(thicknessvar) * 1e4)
    )
    print("Pigment = {} +/- {} wt%".format(w1 * 100, np.sqrt(w1var) * 100))
    print("Kapton = {} +/- {} wt%".format(w2 * 100, np.sqrt(w2var) * 100))

    return energy, muL


if __name__ == "__main__":
    # filename = "/data/id21/inhouse/16dec/hiramFe/spec/16121201.dat"
    # energy1,iodet1,idet1 = gettransmissionxanes(filename,6) # blank
    # energy2,iodet2,idet2 = gettransmissionxanes(filename,9)

    filename = "/data/visitor/hg94/id21/spec/17020201.dat"
    energy1, iodet1, idet1 = gettransmissionxanes(filename, 86)  # blank
    energy2, iodet2, idet2 = gettransmissionxanes(filename, 81)

    b0 = 0
    absorbance_measured = -np.log(idet2 / (iodet2 - b0) * (iodet1 - b0) / idet1)
    energy, absorbance = getcalculated(
        energy2, 2.46, 2.465, 2.49, 2.51, absorbance_measured=absorbance_measured
    )

    fig = plt.figure(1)
    plt.title("{}".format(filename))
    p = plt.plot(energy2, absorbance_measured, label="normalized {}".format(6))
    p = plt.plot(energy, absorbance, label="calculated")
    plt.xlabel("Energy (keV)")
    plt.ylabel("Absorbance")
    plt.legend()
    plt.show()

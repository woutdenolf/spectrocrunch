# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

import pint

ureg = pint.UnitRegistry()


def samplephasor(A, n):
    s = A / np.sqrt(2)
    phasor = np.random.rayleigh(s, n) * np.exp(1j * np.random.uniform(-np.pi, np.pi, n))
    meanabs = s * np.sqrt(np.pi / 2)

    return phasor, meanabs


def approach1():
    wl = ureg.Quantity(1, "angstrom")
    c = ureg.Quantity(1, "c")
    w = 2 * np.pi * (c / wl).to("1/s")
    k0 = ureg.Quantity(np.array([0, 0, 2 * np.pi / wl.to("m").magnitude]), "1/m")
    x = ureg.Quantity(np.array([1, 2, 3]), "mm")
    t = ureg.Quantity(np.linspace(0, (wl / c).to("s").magnitude, 100000), "s")

    # Phasor: real and imaginary part are joint circular Gaussian distributed with zero mean
    # Amplitude: Rayleigh distribution
    # Phase: Uniform in [-pi,pi[
    A = 100

    phasorx, meanabsphasor = samplephasor(A, len(t))
    phasory, meanabsphasor = samplephasor(A, len(t))
    phasor, meanabsphasor = samplephasor(A, len(t))

    # plt.hist(np.abs(phasor), 50)
    # print "Mean amplitude:",meanabsphasor, np.mean(np.abs(phasor))
    # print "Mean phase:",0,np.mean(np.arctan2(np.imag(phasor),np.real(phasor)))
    # plt.show()
    # exit()

    alpha = sum(k0 * x) - w * t

    # Unpolarized
    ux = phasorx * np.exp(1j * alpha)
    uy = phasory * np.exp(1j * alpha)

    plt.figure()
    plt.hist(np.abs(ux) ** 2 + np.abs(uy) ** 2, 50)

    plt.figure()
    plt.plot(np.real(ux), np.real(uy), ".")
    plt.gca().set_aspect("equal", "datalim")

    # Polarized
    beta = np.radians(10.0)
    delta = np.radians(0.0)
    ux = phasor * np.cos(beta) * np.exp(1j * alpha)
    uy = phasor * np.sin(beta) * np.exp(1j * delta) * np.exp(1j * alpha)

    plt.figure()
    plt.hist(np.abs(ux) ** 2, 50)

    plt.figure()
    ux = meanabsphasor * np.cos(beta) * np.exp(1j * alpha)
    uy = meanabsphasor * np.sin(beta) * np.exp(1j * delta) * np.exp(1j * alpha)
    plt.plot(np.real(ux), np.real(uy), linewidth=2)
    chi = np.arctan(np.tan(2 * beta) * np.cos(delta)) / 2.0
    x = np.linspace(-meanabsphasor, meanabsphasor, 2)
    y = np.tan(chi) * x
    plt.plot(x, y, linewidth=2)
    plt.gca().set_aspect("equal", "datalim")  # needed!

    plt.show()


def approach2():

    wl = ureg.Quantity(1, "angstrom")
    c = ureg.Quantity(1, "c")
    w = 2 * np.pi * (c / wl).to("1/s")
    k0 = ureg.Quantity(np.array([0, 0, 2 * np.pi / wl.to("m").magnitude]), "1/m")
    x = ureg.Quantity(np.array([1, 2, 3]), "mm")
    t = ureg.Quantity(np.linspace(0, (wl / c).to("s").magnitude, 100000), "s")

    alpha = sum(k0 * x) - w * t

    A = 100
    phasor, meanabsphasor = samplephasor(A, len(t))

    plt.figure()

    if True:
        P = 1.0
        cosbeta = np.sqrt((1 + P) / 2.0)
        sinbeta = np.sqrt((1 - P) / 2.0)
    else:
        beta = np.radians(20.0)
        cosbeta = np.cos(beta)
        sinbeta = np.sin(beta)
        P = np.cos(2 * beta)

    sin2beta = 2 * sinbeta * cosbeta
    cos2beta = P
    tan2beta = sin2beta / cos2beta

    delta = np.radians(0)

    AA = 1 / (meanabsphasor * np.sin(delta) * cosbeta) ** 2
    CC = 1 / (meanabsphasor * np.sin(delta) * sinbeta) ** 2
    BB = (
        -2
        * np.cos(delta)
        / (meanabsphasor ** 2 * (np.sin(delta)) ** 2 * sinbeta * cosbeta)
    )

    if BB == 0:
        if AA < CC:
            chi1 = 0
        else:
            chi1 = np.pi / 2
    else:
        chi1 = np.arctan((CC - AA - np.sqrt((CC - AA) ** 2 + BB ** 2)) / BB)
    chi2 = np.arctan(tan2beta * np.cos(delta)) / 2.0
    print "Ellipse orientation: ", np.degrees(chi1), np.degrees(chi2)
    chi = chi2

    ux = meanabsphasor * cosbeta * np.exp(1j * alpha)
    uy = meanabsphasor * sinbeta * np.exp(1j * delta) * np.exp(1j * alpha)
    urx = np.real(ux)
    ury = np.real(uy)
    plt.plot(urx, ury, linewidth=2)

    semimajor = meanabsphasor * np.sqrt(
        (1 + np.sqrt(1 - (sin2beta * np.sin(delta)) ** 2)) / 2.0
    )
    semiminor = meanabsphasor * np.sqrt(
        (1 - np.sqrt(1 - (sin2beta * np.sin(delta)) ** 2)) / 2.0
    )

    if delta != 0:
        np.testing.assert_allclose(AA * urx ** 2 + BB * urx * ury + CC * ury ** 2, 1)
        np.testing.assert_allclose(
            semimajor, np.sqrt(2.0 / (AA + CC - np.sqrt((AA - CC) ** 2 + BB ** 2)))
        )
        np.testing.assert_allclose(
            semiminor, np.sqrt(2.0 / (AA + CC + np.sqrt((AA - CC) ** 2 + BB ** 2)))
        )

    p = np.linspace(0, 2 * np.pi, len(t))
    x = semimajor * np.cos(p) * np.cos(chi) - semiminor * np.sin(p) * np.sin(chi)
    y = semimajor * np.cos(p) * np.sin(chi) + semiminor * np.sin(p) * np.cos(chi)
    plt.plot(x, y, linewidth=2)

    x = [0, semimajor * np.cos(chi)]
    y = [0, semimajor * np.sin(chi)]
    plt.plot(x, y, linewidth=2)
    x = [0, semiminor * np.cos(chi + np.pi / 2)]
    y = [0, semiminor * np.sin(chi + np.pi / 2)]
    plt.plot(x, y, linewidth=2)

    plt.gca().set_aspect("equal", "datalim")  # needed!

    plt.show()


if __name__ == "__main__":

    approach2()

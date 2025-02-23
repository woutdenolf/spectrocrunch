import xraylib
import numpy as np
import requests
import re as regex

elements = ["Si", "N"]
n = [3, 4]
density = 3.44  # g/cm^3
energy = 30

wavelength = xraylib.KEV2ANGST / energy * 1e-8  # cm
NA = xraylib.AVOGNUM * 1e24  # atom/mol
re = xraylib.R_E * 1e2  # cm

Z = [xraylib.SymbolToAtomicNumber(e) for e in elements]


def fgen(f, energy):
    def inner(z):
        return f(z, energy)

    return np.vectorize(inner)


fMM = np.vectorize(xraylib.AtomicWeight)
fFii = np.vectorize(fgen(xraylib.Fii))
fmuPE = np.vectorize(fgen(xraylib.CS_Photo_Total))
fmuTot1 = np.vectorize(fgen(xraylib.CS_Total_Kissel))
fmuTot2 = np.vectorize(fgen(xraylib.CS_Total))

nfrac = np.asarray(n, dtype=float)
nfrac /= nfrac.sum()
MM = fMM(Z)  # g/mol
wfrac = nfrac * MM / (nfrac * MM).sum()
N = wfrac * NA / MM * density  # atoms/cm^3


# https://physics.nist.gov/PhysRefData/FFast/Text1995/chap02.html
Fii = fFii(Z)  # Im(f) in e/atom


def fFii(mu, re, wavelength, NA, MM):
    return mu / (-2 * re * wavelength * NA / MM)


print("\nIm(f) = mu/(-2.re.wavelength.NA/MM)")

print(" Z = {}".format(Z))
print(" Im(f)(xraylib) = {}".format(",".join("{:f}".format(x) for x in Fii)))

Fiicalc = fFii(fmuPE(Z))
print(" Im(f)(CS_Photo_Total) = {}".format(",".join("{:f}".format(x) for x in Fiicalc)))

Fiicalc = fFii(fmuTot1(Z))
print(
    " Im(f)(CS_Total_Kissel) = {}".format(",".join("{:f}".format(x) for x in Fiicalc))
)

Fiicalc = fFii(fmuTot2(Z))
print(" Im(f)(CS_Total) = {}".format(",".join("{:f}".format(x) for x in Fiicalc)))


# https://physics.nist.gov/PhysRefData/FFast/Text1995/chap01.html
def fbeta(Imf, re, wavelength, N):
    return -re / (2 * np.pi) * wavelength**2 * sum(N * Imf)


# http://henke.lbl.gov/optical_constants/getdb2.html
url = "http://henke.lbl.gov/cgi-bin/getdb.pl"
params = {
    "Density": density,
    "Formula": "".join("{}{:d}".format(e, _n) for e, _n in zip(elements, n)),
    "Max": energy * 1000,
    "Min": energy * 1000,
    "Npts": 1,
    "Output": "Text File",
    "Scan": "Energy",
}

r = requests.post(url, params)
m = regex.compile('<META.+?CONTENT=".+?URL=(.+?)">')
m = m.search(r.text)

url = "http://henke.lbl.gov{}".format(m.groups()[0])
r = requests.get(url)
r = r.content.split("\n")[2].split(" ")
r = [x for x in r if x]
delta, beta = r[1], r[2]

print("\nn = 1-delta-i.beta")
print("beta = wavelength/(4.pi).density.sum(w.mu)")
print(" beta(Henke) = {}".format(beta))
print(" beta(mu=CS_Photo_Total) = {}".format(fbeta(fFii(fmuPE(Z)))))
print(" beta(mu=CS_Total_Kissel) = {}".format(fbeta(fFii(fmuTot1(Z)))))
print(" beta(mu=CS_Total) = {}".format(fbeta(fFii(fmuTot2(Z)))))

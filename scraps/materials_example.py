# -*- coding: utf-8 -*-

execfile("initcctbx.py")

# Don't use the installed version
import os, sys
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spectrocrunch.materials.compoundfromcif import compoundfromcif
from spectrocrunch.materials.mixture import mixture
from spectrocrunch.materials.types import fraction

import numpy as np
import matplotlib.pyplot as plt
import xraylib

def muplot(mu1,mu2,mixture):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    number = len(mu1)+1
    cmap = plt.get_cmap('jet')
    colors = [cmap(i) for i in np.linspace(0, 1, number)]

    i = 0
    if isinstance(mu1,dict):
        tot = E*0
        for e in mu1:
            tot += mu1[e]
            ax.plot(E, mu1[e], color=colors[i], linestyle='None', marker='o', label='{} ({:.2f} wt%)'.format(e.name,100*mixture[e]))
            i += 1
    else:
        tot = mu1
    ax.plot(E, tot, color=colors[i], linestyle='None', marker='o', label='sum')

    i = 0
    if isinstance(mu2,dict):
        tot = E*0
        for e in mu2:
            tot += mu2[e]
            ax.plot(E, mu2[e], color=colors[i], label='{} ({:.2f} wt%)'.format(e.name,100*mixture[e]))
            i += 1
    else:
        tot = mu2
    ax.plot(E,tot,color=colors[i], label='sum')

    ax.legend(loc='best')
    ax.set_xlabel('Energy (keV)')
    ax.set_ylabel('Mu (cm^2/g)')
    #fig.show() # doesn't block
    plt.show()

if __name__ == '__main__':

    compound1 = compoundfromcif("cif/cinnabar.cif",name='cinnabar')
    compound2 = compoundfromcif("cif/gypsum.cif",name='gypsum')
    mixture = mixture([compound1,compound2],[0.5,0.5],fraction.mass)
    #mixture = compoundfromcif("gypsum.cif",name='gypsum')

    # Shells
    #shells = [xraylib.L1_SHELL,xraylib.L2_SHELL,xraylib.L3_SHELL]
    #shells = [xraylib.K_SHELL]
    shells = [xraylib.K_SHELL]

    # Fluo lines that belong these shells
    fluolines = []
    fluolines += [xraylib.__dict__[s] for s in xraylib.__dict__.keys() if s.endswith("_LINE") and s.startswith("K")]
    #fluolines += [xraylib.__dict__[s] for s in xraylib.__dict__.keys() if s.endswith("_LINE") and s.startswith("L1")]
    #fluolines += [xraylib.__dict__[s] for s in xraylib.__dict__.keys() if s.endswith("_LINE") and s.startswith("L2")]
    #fluolines += [xraylib.__dict__[s] for s in xraylib.__dict__.keys() if s.endswith("_LINE") and s.startswith("L3")]

    # Element
    mixture.markabsorber("S",shells=shells,fluolines=fluolines)

    # Energies
    #steps = np.array([-200, 10, -50, 0.1, 10., 2, 100., 10, 800.])
    #steps = np.array([-100, 10, -50, 1, 50., 2, 100., 10, 150.])    
    #E = mixture.get_energy(steps)
    E = np.linspace(2.46,2.52,120)

    # Cross-sections
    #mu1 = mixture.mass_att_coeff(E,decomposed=False,fine=True)
    #mu2 = mixture.mass_att_coeff(E,decomposed=False)
    mu1 = mixture.partial_mass_abs_coeff(E,decomposed=False,fine=True,refresh=True)
    mu2 = mixture.partial_mass_abs_coeff(E,decomposed=False)

    # Plot
    muplot(mu1,mu1,mixture)

    #tmp = raw_input("Press enter to continue ...")


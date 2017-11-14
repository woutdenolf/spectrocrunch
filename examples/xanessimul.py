
import os, sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spectrocrunch.materials.compoundfromcif import compoundfromcif

import matplotlib.pyplot as plt
import numpy as np
import xraylib

if __name__ == '__main__':
    c = compoundfromcif("cif/calcite.cif",name='calcite')
    
    energy = np.linspace(3.9,4.4,100)
    c.markabsorber('Ca',xraylib.K_SHELL,fluolines=[])
    mu1 = c.mass_att_coeff(energy,fine=False)
    mu2 = c.mass_att_coeff(energy,fine=True,refresh=False)
    
    plt.plot(energy,mu1*c.density,'o')
    plt.plot(energy,mu2*c.density,'o')
    plt.show()
    

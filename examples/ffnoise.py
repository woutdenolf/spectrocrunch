import os, sys
sys.path.insert(1,'/data/id21/inhouse/wout/dev/SpectroCrunch')

from spectrocrunch.io.spec import spec
import matplotlib.pyplot as plt
import numpy as np

#from spectrocrunch.materials.compoundfromformula import compoundfromformula as compound
from spectrocrunch.materials.compoundfromcif import compoundfromcif as compound
from spectrocrunch.materials.mixture import mixture
from spectrocrunch.materials.types import fractionType

import xraylib
import silx.math.fit as fit

def gettransmissionxanes(filename,scannumber):
    f = spec(filename)
    result = f.getdata2(scannumber, ["arr_energyM","arr_iodet","arr_idet"])
    energy = result[:,0]
    iodet = result[:,1]
    idet = result[:,2]
    return energy,iodet,idet

def muL_combination(x,thickness,w1):
    energy, b1, b2, rho1, rho2 = x
    w2 = 1-w1
    rho = 1/(w1/rho1+w2/rho2)
    return thickness*rho*(w1*b1+w2*b2)

def getcalculated(energy,e1,e2,absorbance_measured=None):
    #compound1 = compound("CaCO3",2.71,name='calcite')
    #compound2 = compound("Fe2O3",5.3,name='hematite')
    compound1 = compound("calcite",name='calcite')
    compound2 = compound("hematite",name='hematite')
    m = mixture([compound1,compound2],[0.5,0.5],fractionType.weight)
    print(m.elemental_weightfractions())
    print(m.elemental_molefractions())

    mask = (energy <= e1) | (energy >= e2)

    mu = m.mass_att_coeff(energy[mask],decomposed=True)
    b1 = mu['calcite']*2
    b2 = mu['hematite']*2

    thickness = 5e-4
    w1 = 0.5
    w2 = 0.5

    if absorbance_measured is not None:
        x = energy[mask], b1, b2, compound1.density, compound2.density
        y = absorbance_measured[mask]
        ysigma = np.sqrt(y)

        p0 = (thickness,w1)
        constraints = [[fit.CFREE,0,0],[fit.CQUOTED,0,1]]

        p, cov_matrix, info = fit.leastsq(muL_combination, x, y, p0, sigma=ysigma, constraints=constraints, full_output=True)
        var = np.diag(cov_matrix)

        S = np.diag(1/np.sqrt(var))
        cor_matrix = S.dot(cov_matrix).dot(S)
        print("R(d,w1) = {}".format(cor_matrix[0,1]))
        

        thickness,w1 = p
        thicknessvar,w1var = var
        w2 = 1-w1
        w2var = w1var

    m = mixture([compound1,compound2],[w1,w2],fractionType.weight)
    muL = m.mass_att_coeff(energy,decomposed=False)*m.density()*thickness

    print("Density = {} g/cm^3".format(m.density()))
    print("Thickness = {} +/- {} um".format(thickness*1e4,np.sqrt(thicknessvar)*1e4))
    print("Calcite = {} +/- {} wt%".format(w1*100,np.sqrt(w1var)*100))
    print("Hematite = {} +/- {} wt%".format(w2*100,np.sqrt(w2var)*100))
    
    return energy,muL
    

if __name__ == '__main__':
    filename = "/data/id21/inhouse/16dec/hiramFe/spec/16121201.dat"
    energy1,iodet1,idet1 = gettransmissionxanes(filename,6)
    energy2,iodet2,idet2 = gettransmissionxanes(filename,9)

    b0 = 300
    absorbance_measured = -np.log(idet2/(iodet2-b0)*(iodet1-b0)/idet1)
    energy,absorbance = getcalculated(energy2,7.1,7.25,absorbance_measured=absorbance_measured)

    fig = plt.figure(1)
    plt.title("{}".format(filename))
    p = plt.plot(energy2,absorbance_measured,label="normalized {}".format(6))
    p = plt.plot(energy,absorbance,label="calculated")
    plt.xlabel("Energy (keV)")
    plt.ylabel("Absorbance")
    plt.legend()
    plt.show()



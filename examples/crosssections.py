
import os, sys
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spectrocrunch.materials.compoundfromformula import compoundfromformula as compound

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

def plot(c,path=None):
    energy = np.linspace(2,30,500)

    mu = c.mass_att_coeff(energy)
    muabs = c.mass_abs_coeff(energy)
    mucompton = c.compton_cross_section(energy)
    murayleigh = c.rayleigh_cross_section(energy)

    thickness = 50
    Ptrans = np.exp(-mu*c.density*thickness*1e-4)
    Patt = 1-Ptrans
    Pabs = Patt * muabs/mu
    Pcompton = Patt * mucompton/mu
    Prayleigh = Patt * murayleigh/mu

    plt.figure(1)
    plt.clf()
    plt.plot(energy,Ptrans*100,label='No interaction')
    plt.plot(energy,Pabs*100,label='Photo-ionization')
    plt.plot(energy,Pcompton*100,label='Inellastic scattering')
    plt.plot(energy,Prayleigh*100,label='Elastic scattering')
    #plt.plot(energy,muabs/mu*100,label='Photo-ionization')
    #plt.plot(energy,mucompton/mu*100,label='Inellastic scattering')
    #plt.plot(energy,murayleigh/mu*100,label='Elastic scattering')

    plt.plot()
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_ylim([1e-3,1e2])
    ax.set_xlabel('Energy (keV)')
    ax.set_ylabel('Probability @ 50 $\mu$m (%)')
    #ax.set_ylabel('Probability (%)')
    plt.legend(loc=4)
    

    if path is not None:
        plt.savefig(os.path.join(path,"{}.cs.png".format(c)))

    plt.figure(2)
    plt.clf()
    thickness = np.linspace(10,50,9)
    ttheta = None
    colors = iter(cm.gist_rainbow(np.linspace(0, 1, len(thickness))))

    for t in thickness:
        scatyield = c.density*murayleigh*np.exp(-t*1e-4*c.density*mu)
        
        if ttheta is None:
            scatyield *= t*1e-4
        else:
            a = mu*c.density*(1-1/np.cos(ttheta*np.pi/180))
            print (1-np.exp(-a*t*1e-4))/a
            print t*1e-4
            scatyield *= (1-np.exp(-a*t*1e-4))/a

        plt.plot(energy,scatyield*100,label="{} $\mu$m".format(t),color=next(colors))

    ax = plt.gca()
    ax.set_xlabel('Energy (keV)')
    ax.set_ylabel('Forward scattering yield (%)')
    plt.legend(loc=1, fontsize = 'x-small')

    if path is not None:
        plt.savefig(os.path.join(path,"{}.syield.png".format(c)))

    if path is None:
        plt.show()

if __name__ == '__main__':
    
    cmps = [compound("Fe2O3",5.26),\
            compound("CaCO3",2.71),\
            compound("Pb3(CO3)2(OH)2",6.8)]

    for c in cmps:
        plot(c,path="/data/id21/inhouse/wout/tmp/cs")


    
    


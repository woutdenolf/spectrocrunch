
import sys
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spectrocrunch.io.edf import get2Dscancoordinates as getcoord

import fabio
import os
import numpy as np
import matplotlib.pyplot as plt

def getpos(path,sample,scannr,zapctr):
    f = os.path.join(path,sample,"zap","{}_arr_{}_{:04}_0000.edf".format(sample,zapctr,scannr))
    f = fabio.open(f)
    return getcoord(f.header,scanlabel="Title")

def getencoder(path,sample,scannr,mot):
    f = os.path.join(path,sample,"zap","{}_arr_{}_{:04}_0000.edf".format(sample,mot,scannr))
    f = fabio.open(f)
    return f.data

def stepspermotunit(pos,epos,n):
    d = np.median(pos[1:]-pos[0:-1])
    de = np.asarray([np.median(epos[i,1:]-epos[i,0:-1]) for i in range(n)])
    return np.median(de/d)

if __name__ == '__main__':
    path = "/data/id21/inhouse/17mar/Hercules"
    sample = "As1test"
    scannr = 1

    pos = getpos(path,sample,scannr,"iodet")
    efast = getencoder(path,sample,scannr,pos[0]["name"])
    eslow = getencoder(path,sample,scannr,pos[1]["name"])
    pfast = pos[0]["data"]
    pslow = pos[1]["data"]
    nfast = len(pos[0]["data"])
    nslow = len(pos[1]["data"])
    fastname = pos[0]["name"]
    slowname = pos[1]["name"]

    # Fast dimension horizontal
    if pos[0]["name"]=="samz":
        efast = efast.T
        eslow = eslow.T
   
    # Motor resolution
    resfast = stepspermotunit(pfast,efast,nslow)
    resslow = stepspermotunit(pslow,eslow.T,nfast)
    print("{}: {} steps/micron".format(fastname,resfast/1000))
    print("{}: {} steps/micron".format(slowname,resslow/1000))


    plt.figure(1)
    for i in range(nslow):
        plt.plot(pfast,efast[i,:])

    plt.figure(2)
    for i in range(nfast):
        plt.plot(pslow,eslow[:,i])

    plt.show()



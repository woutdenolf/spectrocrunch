# -*- coding: utf-8 -*-

filename = "/data/id16b/inhouse1/comm_17jan/restart/sofc/26jan/6100h_fluoXAS_0/results/6100h_fluoXAS_0/test.h5"
name = "/detectorsum/Ni-K"
new = "/detectorsum/Ni-K_norm"

specfile = "/data/id16b/inhouse1/comm_17jan/ma3257/align.spec"
specscannumber = 33
energyname = "energy"
fluxname = "flux_It"







###### Import libraries ######
from spectrocrunch.io.spec import spec
from spectrocrunch.h5stacks.get_hdf5_imagestacks import get_hdf5_imagestacks 
from spectrocrunch.math.interpolate import extrap1d
from scipy.interpolate import interp1d
import spectrocrunch.io.nexus as nexus
import h5py
import numpy as np
import matplotlib.pyplot as plt

###### Get flux from spec file ######
fspec = spec(specfile)
data,info = fspec.getdata(specscannumber,[energyname,fluxname])
energy = data[:,0]
Inorm = data[:,1]
Inorm[:] = 2.
fluxfunc = extrap1d(interp1d(energy,Inorm))

###### Get stack to normalize ######
stacks, axes = get_hdf5_imagestacks(filename,["detectorsum"])

###### Add normalized dataset ######
fh5 = h5py.File(filename)
stackdim = 2

stackenergy = fh5[axes[stackdim]["fullname"]][:]
stacknorm = fluxfunc(stackenergy)

plt.plot(energy,Inorm,"-",label="Spec")
plt.plot(stackenergy,stacknorm,"or",label="fluoXAS")
plt.xlabel("Energy (keV)")
plt.ylabel("Flux")
plt.title("Flux used for normalization")
plt.legend()
plt.show()

tmp = [s for s in new.split("/") if len(s)>0]
nxdatagrp = nexus.newNXdata(fh5[tmp[0]],'/'.join(tmp[1:]),"")
dset = nexus.createNXdataSignal(nxdatagrp,shape=fh5[name]["data"][:].shape,chunks = True,dtype = fh5[name]["data"][:].dtype)
nexus.linkaxes(fh5,axes,[nxdatagrp])

data = fh5[name]["data"][:]
s = np.array(data.shape)
#stackdim = np.where(s==stackenergy.size)[0][0]

s = np.delete(s,stackdim)
if stackdim == 0:
    snew = (s[1],s[0],stacknorm.size)
elif stackdim == 1:
    snew = (s[1],stacknorm.size,s[0])
else:
    snew = (stacknorm.size,s[1],s[0])

data /= np.tile(stacknorm,(s[0]*s[1],1)).T.reshape(snew).T
dset[:] = data

fh5.close()



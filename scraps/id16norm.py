# -*- coding: utf-8 -*-

filename = "/data/id16b/inhouse1/comm_17jan/restart/sofc/26jan/6100h_fluoXAS_0/results/6100h_fluoXAS_0/test.h5"
name = ["/detectorsum/Ni-K_norm", "/counters/arr_srcurr"]
new = "/detectorsum/Ni-K_norm2"


###### Import libraries ######
import spectrocrunch.io.nexus as nexus
from spectrocrunch.h5stacks.get_hdf5_imagestacks import get_hdf5_imagestacks
import h5py
import numpy as np

###### Normalize one dataset by another ######
stacks, axes = get_hdf5_imagestacks(filename, ["detectorsum"])

f = h5py.File(filename)
tmp = [s for s in new.split("/") if len(s) > 0]
nxdatagrp = nexus.newNXdata(f[tmp[0]], "/".join(tmp[1:]), "")
dset = nexus.createNXdataSignal(
    nxdatagrp,
    shape=f[name[0]]["data"][:].shape,
    chunks=True,
    dtype=f[name[0]]["data"][:].dtype,
)
nexus.linkaxes(f, axes, [nxdatagrp])

dset[:] = f[name[0]]["data"][:] / f[name[1]]["data"][:]


f.close()

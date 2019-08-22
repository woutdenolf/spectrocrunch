# -*- coding: utf-8 -*-

import os,sys
sys.path.insert(1,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spectrocrunch.xrf.parse_xia import parse_xia_esrf
from spectrocrunch.io import xiaedf

from time import time

import shutil
import numpy as np

def expand(data,ndet=6):
    data = data[np.newaxis,np.newaxis,:,np.newaxis]
    data = np.repeat(data,40,axis=0)
    data = np.repeat(data,43,axis=1)
    data = np.repeat(data,ndet,axis=3)
    return data

def gentestdata(datadir,scanname,scannumber):
    data = expand(np.arange(4096,dtype=np.int32))
    stats = expand(np.asarray([0,0,2,1,0,0],dtype=np.int32))
    o = xiaedf.xiaimage_number(datadir,scanname,scannumber)
    xialabels = ["xia{:02d}".format(i) for i in range(data.shape[-1])]
    o.save(data,xialabels,stats=stats)

if __name__ == '__main__':
    test = True

    if test:
        datadir = "/data/id21/inhouse/wout/laurence/test"
        scanname = "test"
        scannumber = 1
        if os.path.isdir(datadir):
            shutil.rmtree(datadir)
        gentestdata(datadir,scanname,scannumber)
    else:
        datadir = "/data/id21/inhouse/wout/laurence/Bloc3L_SynL_cell5_fine2"
        scanname = "Bloc3L_SynL_cell5_fine2"
        scannumber = 1

    gendata = True
    deadtime = True
    add = True
    norm = False
    exclude_detectors=[]
    outdir = "/data/id21/inhouse/wout/laurence/cor/"

    outdir1 = os.path.join(outdir,"1")
    outdir2 = os.path.join(outdir,"2")
    outname = "cor"
    normctr = "zap_it"

    if gendata:
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)

        # Normal parsing
        t0 = time()
        parse_xia_esrf(datadir,scanname,scannumber,outdir1,outname,exclude_detectors=exclude_detectors,\
                        deadtime=deadtime,add=add)
        print time()-t0

        outname = "{}_{}".format(scanname,outname)

        # Advanced parsing
        t0 = time()
        mapin = xiaedf.xiaimage_number(datadir,scanname,scannumber)
        mapin.exclude_detectors(["xiaS0"]+["xia{:02d}".format(i) for i in exclude_detectors])
        mapin.dtcor(deadtime)
        if norm:
            mapin.norm(normctr)
        else:
            mapin.norm(None)
        mapin.detectorsum(add)
        mapout = xiaedf.xiaimage_number(outdir2,outname,scannumber)
        data = mapin.data
        
        if add:
            xialabels = ["xiaS1"]
        else:
            xialabels = ["xia{:02d}".format(i) for i in range(data.shape[-1])]
        mapout.save(data,xialabels)
        print time()-t0
        
    
    map1 = xiaedf.xiaimage_number(outdir1,outname,scannumber)
    map2 = xiaedf.xiaimage_number(outdir2,outname,scannumber)
    data1 = map1.data
    data2 = map2.data

    if test:
        m = 1
        if deadtime:
            m *= 2
        if add:
            m *= 6
            ndet = 1
        else:
            ndet = 6

        data3 = expand(np.arange(4096,dtype=np.int32)*m,ndet=ndet)
        np.testing.assert_allclose(data1,data3)
        np.testing.assert_allclose(data2,data3)
    else:
        np.testing.assert_allclose(data1,data2)


    


    

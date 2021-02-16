# -*- coding: utf-8 -*-

import os, sys

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spectrocrunch.align.alignElastix import alignElastix as align

# from spectrocrunch.align.alignSimple import alignMax as align
# from spectrocrunch.align.alignSift import alignSift as align

import h5py

if __name__ == "__main__":
    path = "/data/visitor/ls2497/id21/XANES from Diamond/"
    source = os.path.join(path, "Xanes1Co.h5")
    dest = os.path.join(path, "Xanes1Co.aligned.h5")
    stackdim = 2
    sourcelist = ["/exchange/data"]
    destlist = ["/exchange/data"]
    roi = ((0, -1), (5, 50))

    a = 15
    b = 58
    c = 20
    hin = h5py.File(source, "r")
    source = [hin["/exchange/data"][..., a:b]]
    sourcelist = "None"

    os.remove(dest)
    o = align(
        source,
        sourcelist,
        dest,
        destlist,
        "",
        stackdim=stackdim,
        overwrite=True,
        plot=True,
    )
    o.align(0, refimageindex=c, onraw=True, pad=False, crop=True, roi=roi)

    hout = h5py.File(dest, "a")
    hout["/exchange/energy"] = hin["/exchange/energy"][a:b]
    hin.close()
    hout.close()

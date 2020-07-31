# -*- coding: utf-8 -*-

import os, sys

sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from spectrocrunch.process.id21_xas import processEnergySynchronized as process

if __name__ == "__main__":
    path = os.path.dirname(os.path.abspath(__file__))

    specfile = "/data/visitor/hg79/id21/spec/16050901.dat"
    specnumbers = [[268, 269, 270], [272, 273, 274]]

    detectorcfg = os.path.join(path, "testdata", "xrfxanes", "id21", "hg79.cfg")

    destpath = os.path.join(path, "testresults", "hg79")
    destradix = "test"
    process(
        specfile,
        specnumbers,
        destpath,
        destradix,
        detectorcfg,
        showelement="Pb M",
        sourcedir=None,
        dtcor=True,
        fastfitting=True,
    )

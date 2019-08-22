# -*- coding: utf-8 -*-

from spectrocrunch.process.id21_fluoxas import process

if __name__ == '__main__':

    # Data
    scanname = "fXANES5"
    scannumbers = range(1,173) # from a to b-1 !!!!!
    sourcepath =  "/data/id21/inhouse/15apr/Hiram/fXANES5"

    # Fitting
    cfgfile = "/data/id21/inhouse/15apr/Hiram/results/fXANES5.cfg" # or None

    # Results
    destpath = "/data/id21/inhouse/15apr/Hiram/results/fXANES5"

    skippreprocessing = True

    # Alignment
    alignmethod = "fft" #None, fft, sift, elastix
    alignreference = "Ca-K"

    # Do processing
    process(sourcepath,destpath,scanname,scannumbers,cfgfile,alignmethod,alignreference,skippre=skippreprocessing)


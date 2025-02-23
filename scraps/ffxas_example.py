import os
import logging

from spectrocrunch.process.id21_ffxas import process

logging.getLogger("spectrocrunch").setLevel(logging.INFO)

if __name__ == "__main__":

    path = os.path.dirname(os.path.abspath(__file__))

    # Raw data
    sourcepath = "/data/id21/inhouse/wout/dev/SpectroCrunch/examples/testdata/ff/id21/7913_50RH_ffCd"
    radix = "map1"
    rebin = (1, 1)
    skippreprocessing = True

    # Result
    destpath = os.path.join(path, "testresults", "7913_50RH_ffCd")

    # Normalization
    skipnormalization = True

    # Alignment
    alignmethod = "elastix"
    refimageindex = 0
    crop = True
    roi = ((100, 260), (160, 350))
    plot = True

    # Process
    process(
        sourcepath,
        destpath,
        radix,
        "",
        rebin,
        alignmethod,
        skippre=skippreprocessing,
        skipnormalization=skipnormalization,
        refimageindex=refimageindex,
        crop=crop,
        roialign=roi,
        plot=plot,
    )

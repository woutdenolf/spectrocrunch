import os
import sys
import logging

from spectrocrunch.process.id21_fluoxas import process

logging.getLogger("spectrocrunch").setLevel(logging.INFO)

if __name__ == "__main__":
    # defaultstack('/mntdirect/_data_id21_inhouse/wout/dev/SpectroCrunch/examples/testresults/fXANES5/fXANES5.norm.h5',['/detector0/Ca-K'],'Ca-K')
    # exit()

    path = os.path.dirname(os.path.abspath(__file__))

    example = "example6"

    cfgfile = None
    skippreprocessing = False
    skipnormalization = False
    skipalign = False
    dtcor = False
    default = None
    crop = True
    roialign = None
    plot = True
    bff = False

    if example == "example1":
        scanname = ["fXANES5"]
        scannumbers = [[1, 50, 100, 170]]  # [range(1,173)]
        sourcepath = [os.path.join(path, "testdata", "xrfxanes", "id21", scanname[0])]
        destpath = os.path.join(path, "testresults", scanname[0])
        cfgfile = os.path.join(
            path, "testdata", "xrfxanes", "id21", scanname[0] + ".cfg"
        )

        skippreprocessing = False
        skipnormalization = False

        alignmethod = "fft"  # None, fft, sift, elastix
        alignreference = "Ca-K"
        refimageindex = 0  # None for pair-wise alignment
        crop = True
        roialign = ((10, 60), (10, 60))

        default = "Ca-K"

    elif example == "example2":
        scanname = ["saphir"]
        scannumbers = [range(1, 143)]
        sourcepath = [os.path.join(path, "testdata", "xrfxanes", "id21", scanname[0])]
        destpath = os.path.join(path, "testresults", scanname[0])
        cfgfile = os.path.join(
            path, "testdata", "xrfxanes", "id21", scanname[0] + ".cfg"
        )
        scanname[0] += "1_xrf_ff"

        skippreprocessing = True
        skipnormalization = False

        alignmethod = None  # None, fft, sift, elastix
        alignreference = "xmap_x3c"
        refimageindex = 0  # None for pair-wise alignment

    elif example == "example3":
        scanname = ["roi4_p1", "roi4_p1"]
        scannumbers = [range(1, 54), range(1, 63 - 10)]
        sourcepath = [
            "/data/visitor/in946/id21/fXAS_sample1/roi4_p1_fluoXAS_7",
            "/data/visitor/in946/id21/fXAS_sample1/roi4_p1_fluoXAS_8",
        ]
        destpath = os.path.join(path, "testresults", "roi4_p1")
        cfgfile = None

        skippreprocessing = False
        skipnormalization = False

        alignmethod = "elastix"  # None, fft, sift, elastix
        alignreference = "xmap_x3c"
        refimageindex = 6  # None for pair-wise alignment

    elif example == "example4":
        scanname = ["fXANESpl"]
        scannumbers = [
            range(1, 5) + range(50, 53)
        ]  # [range(1,84)+range(85,130)] # from a to b-1 !!!!!
        sourcepath = ["/data/id21/inhouse/15fev/Hiram/fXANESplasma"]

        # Fitting
        cfgfile = None  # "/data/id21/inhouse/15apr/Hiram/results/fXANES5.cfg" # or None

        # Results
        destpath = os.path.join(path, "testresults", scanname[0])

        skippreprocessing = False
        skipnormalization = False

        # Alignment
        alignmethod = "max"  # None, "fft", "sift", "elastix"
        alignreference = "xmap_x3c"
        refimageindex = 1  # None or a number

    elif example == "example5":
        scanname = ["KMS"]
        scannumbers = [[5, 6, 7, 8]]
        sourcepath = [os.path.join(path, "testdata", "xrfxanes", "id16b", "KMSa")]

        cfgfile = [
            "KMSa_En4_Det2.cfg",
            "KMSa_En4_Det3.cfg",
            "KMSa_En4_Det4.cfg",
            "KMSa_En4_Det6.cfg",
            "KMSa_En4_Det7.cfg",
        ]
        cfgfile = [
            os.path.join(path, "testdata", "xrfxanes", "id16b", s) for s in cfgfile
        ]

        destpath = os.path.join(path, "testresults", "KMSa")

        dtcor = True

        skippreprocessing = True
        skipnormalization = False

        alignmethod = "fft"  # None, fft, sift, elastix
        alignreference = "Fe-K"
        refimageindex = 0  # None for pair-wise alignment
        crop = True

        default = "Fe-K"

    elif example == "example6":
        bff = True

        sourcepath = [os.path.join(path, "testdata", "ff", "id21", "7913_50RH_ffCd")]
        scanname = ["map1"]
        destpath = os.path.join(path, "testresults", "7913_50RH_ffCd ")

        skippreprocessing = False
        skipnormalization = False
        skipalign = False

        alignmethod = "fft"  # None, fft, sift, elastix
        refimageindex = 0  # None for pair-wise alignment
        crop = True

        rebin = (1, 1)
        roiraw = ((300, 650), (900, 1100))
        roialign = ((1, -2), (1, -2))
        roiresult = ((10, -20), (11, -21))

    else:
        sys.exit()

    if bff:
        process(
            sourcepath,
            destpath,
            scanname,
            "",
            rebin,
            alignmethod,
            refimageindex=refimageindex,
            skippre=skippreprocessing,
            skipnormalization=skipnormalization,
            roialign=roialign,
            crop=crop,
            plot=plot,
            skipalign=skipalign,
            roiresult=roiresult,
            roiraw=roiraw,
        )
    else:
        process(
            sourcepath,
            destpath,
            scanname,
            scannumbers,
            cfgfile,
            alignmethod,
            alignreference,
            default=default,
            refimageindex=refimageindex,
            skippre=skippreprocessing,
            skipnormalization=skipnormalization,
            roialign=roialign,
            crop=crop,
            plot=plot,
            dtcor=dtcor,
        )

import os
import h5py
from spectrocrunch.io import xiaedf
from spectrocrunch.io import edf
import matplotlib.pyplot as plt
from PyMca5.PyMcaIO.ConfigDict import ConfigDict


def checkdtcor():
    xiamapcor = xiaedf.xiaimage_number(
        "/data/visitor/es934/id16b/standard/srm1832/DT_corr", "srm1832_corr", 1
    )
    xiamapcor.include_detectors = [4]
    xiamapcor.dtcor(False)
    d1 = xiamapcor.data

    if False:
        xiamaporg = xiaedf.xiaimage_number(
            "/data/visitor/es934/id16b/standard/srm1832", "srm1832", 1
        )
        xiamaporg.dtcor(True)
    else:
        xiamaporg = xiaedf.xiaimage_number(
            "../tmpresults/srm1832/srm1832.map1_pymca.3/xrfspectra/", "srm1832_dtcor", 0
        )
        xiamapcor.dtcor(False)
    xiamaporg.include_detectors = [4]
    d2 = xiamaporg.data

    spectrum1 = d1.sum(axis=(0, 1, 3))
    spectrum2 = d2.sum(axis=(0, 1, 3))
    plt.plot(spectrum1)
    plt.plot(spectrum2)
    plt.yscale("log")
    plt.figure()
    plt.plot(spectrum1 / spectrum2)
    plt.show()


def comparedict(d1, d2, nodes=tuple()):
    if isinstance(d1, dict):
        assert isinstance(d2, dict), nodes
        keys1 = set(d1.keys())
        keys2 = set(d2.keys())
        if keys1 != keys2:
            print(f"{nodes}: missing={keys1-keys2}, unexpected={keys2-keys1}")
        for k in d1:
            if k in d2:
                comparedict(d1[k], d2[k], nodes=nodes + (k,))
    else:
        if d1 != d2:
            print(f"{nodes}: expected={d1}, actual={d2}")


def checkresult():

    # pymca gui DTcor:
    img1 = "/data/visitor/es934/id16b/standard/srm1832/DT_corr/IMAGES/srm1832_corr_xia04_0001_0000_0000_to_0050_V_K.edf"
    fcfg1 = "/data/visitor/es934/id16b/standard/srm1832/DT_corr/config_lagaffe_elt_04_batch.cfg"

    # pymca gui not DTcor, no weights:
    img1 = "../IMAGES_NO_DT/images_V_K.edf"
    fcfg1 = "../IMAGES_NO_DT/images.cfg"

    # pymca gui DTcor, no weights:
    img1 = "../IMAGES/images_V_K.edf"
    fcfg1 = "../IMAGES/images.cfg"

    # SPECTROCRUNCH
    # proc = "pymca.1" # NO DT
    proc = "pymca.2"  # DT
    # proc = "pymca.3" # DT ON SPECTRA

    # img2 = f"../tmpresults/srm1832/srm1832.map1_{proc}/pymcaresults/srm1832_xia04_0001_0000_V_K.edf"
    img2 = f"../tmpresults/srm1832/edfresults/srm1832.map1/{proc}/xia04_V_K.edf"

    fcfg2 = f"../tmpresults/srm1832/srm1832.map1_{proc}/pymcaresults/"
    if proc == "pymca.3":
        fcfg2 += "srm1832_dtcor_xia04_0000_0000.cfg"
    else:
        fcfg2 += "srm1832_xia04_0001_0000.cfg"

    img1 = edf.edfimage(img1).data
    img2 = edf.edfimage(img2).data
    cfg1 = ConfigDict()
    cfg1.read(os.path.abspath(fcfg1))
    cfg2 = ConfigDict()
    cfg2.read(os.path.abspath(fcfg2))

    comparedict(cfg1, cfg2)

    xiamaporg = xiaedf.xiaimage_number(
        "/data/visitor/es934/id16b/standard/srm1832", "srm1832", 1
    )
    xiamaporg.include_detectors = [4]
    xiamaporg.onlyicrocr(True)
    stats = xiamaporg.stats
    dtcor = xiaedf.deadtimecorrector(stats[..., 0, 0], stats[..., 1, 0])

    plt.figure()
    im = plt.imshow(img1 - img2)
    plt.colorbar(im)


checkresult()
# checkdtcor()
plt.show()

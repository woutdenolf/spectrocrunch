from spectrocrunch.visualization.id21_scanoverlap import plot

if __name__ == "__main__":

    specfilename = "/data/id21/store/backup_visitor/2016/hg79/id21/spec/16051103.dat"
    specnumbers = range(5, 110, 4)
    specnumbers = [
        5,
        9,
        13,
        17,
        21,
        25,
        29,
        33,
        37,
        41,
        45,
        49,
        53,
        57,
        61,
        65,
        69,
        73,
        77,
        81,
        85,
        89,
        93,
        97,
        101,
        105,
        109,
    ]

    hdf5filename = (
        "/data/id21/inhouse/xrfalign/hg79/virginsphere/monnight_maps.align.h5"
    )
    offsamz = 0
    offsamy = 0

    grps = {}
    i = 0
    grps[i] = {"path": "/detector0/Pb-M", "ind": 3, "lo": 0.05, "hi": 0.00}
    i += 1
    grps[i] = {"path": "/detector0/Cl-K", "ind": 3, "lo": 0.05, "hi": 0.50}
    i += 1
    grps[i] = {"path": "/detector0/S-K", "ind": 3, "lo": 0.05, "hi": 0.50}
    i += 1

    plot(
        hdf5filename,
        grps,
        specfilename,
        specnumbers,
        offsamy,
        offsamz,
        transpose=False,
        flipvert=True,
        fliphor=False,
        defaultorigin=False,
        showlabels=False,
        color="#ffffff",
        printpos=False,
    )

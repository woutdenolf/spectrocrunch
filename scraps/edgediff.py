import os
import glob
import re
import numpy as np
import fabio

from spectrocrunch.align.alignSift import alignSift as align


def edgeimagepairs(path, radix, threshold=0):
    """path/radix_[0-9]_edge_[0-9].[0-9].edf"""
    base = os.path.join(path, radix)
    pattern = re.compile(base + r"_([0-9]+)_edge_([0-9]+\.?[0-9]*).edf")

    files = sorted(glob.glob(base + "*.edf"))
    files = [f for f in files if pattern.match(f)]

    edgefiles = []
    numbers = []
    energies = []
    for f in files:
        mobj = pattern.match(f)
        if mobj is not None:
            # nr = int(mobj.group(1))
            nr = os.path.getctime(f)
            if nr >= threshold:
                numbers.append(nr)
                energies.append(float(mobj.group(2)))
                edgefiles.append(os.path.basename(f))

    edgefiles = [f for (num, en, f) in sorted(zip(numbers, energies, edgefiles))]

    if len(edgefiles) % 2 != 0:
        edgefiles = edgefiles[:-1]
        numbers = numbers[:-1]
        energies = energies[:-1]

    return edgefiles, energies


def parsepairs(inpath, outpath, edgefiles, energies):

    outputstack = [np.zeros(1, dtype=np.float32)]

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    order = 0

    for i in range(0, len(edgefiles), 2):
        order += 1
        if energies[i] > energies[i + 1]:
            a = i + 1
            b = i
        else:
            a = i
            b = i + 1

        # Prepare alignment: 1 stack of 2 images
        print("Align {} and {}...".format(edgefiles[a], edgefiles[b]))
        o = align(
            inpath,
            [[edgefiles[a], edgefiles[b]]],
            outputstack,
            None,
            None,
            stackdim=2,
            overwrite=True,
        )

        # Align
        o.align(0, onraw=True, extend=False)

        # Save output
        outfile = str.split(edgefiles[a], "_")
        outfile = "_".join(outfile[0:2])

        outfilea = os.path.join(
            outpath, "{:04d}_{}_{}.edf".format(order, "before", outfile)
        )
        outfileb = os.path.join(
            outpath, "{:04d}_{}_{}.edf".format(order, "after", outfile)
        )
        outfilec = os.path.join(
            outpath, "{:04d}_{}_{}.edf".format(order, "diff", outfile)
        )

        ediff = outputstack[0][:, :, 1] - outputstack[0][:, :, 0]

        o = fabio.edfimage.edfimage(data=outputstack[0][:, :, 0])
        o.write(outfilea)
        o = fabio.edfimage.edfimage(data=outputstack[0][:, :, 1])
        o.write(outfileb)
        o = fabio.edfimage.edfimage(data=ediff)
        o.write(outfilec)


if __name__ == "__main__":
    path = "/data/id21/inhouse/16apr/ihch1067/papysample9best/ff/images"
    outpath = "/data/id21/inhouse/16apr/ihch1067/ffproc/edgediff/papysample9best"
    threshold = 0  # os.path.getctime(os.path.join(path,"ff_0001_edge_7.120.edf"))-1e-7

    edgefiles, energies = edgeimagepairs(path, "ff", threshold=threshold)

    parsepairs(path, outpath, edgefiles, energies)

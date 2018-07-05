# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from ..align.alignElastix import alignElastix
from ..align.alignSift import alignSift
from ..align.alignFFT import alignFFT
from ..align.alignSimple import alignMax
from ..align.alignSimple import alignMin
from ..align.alignSimple import alignCentroid
from ..io import nexus

import h5py

def align_hdf5_imagestacks(filein,stacks,axes,stackdim,fileout,alignmethod,refdataset,\
            refimageindex=None,overwrite=False,info=None,copygroups=None,crop=False,roi=None,plot=True):
    """
    Args:
        filein (str)
        stacks (list(str)): list of nxdata paths
        axes (dict)
        stackdim (int)
        fileout (str)
        alignmethod (str)
        refdataset (str)
        refimageindex (Optional(int))
        overwrite (Optional(bool))
        info (dict)
        copygroups (list(str)): hdf5 dataset paths
        crop (Optional(bool)): crop aligned images, padded otherwise
        roi(Optional(array-like)): use ony part of the image to align
        plot (Optional(bool))
    """

    # Alignment method
    if alignmethod=="sift":
        alignclass = alignSift
    elif alignmethod=="elastix":
        alignclass = alignElastix
    elif alignmethod=="fft":
        alignclass = alignFFT
    elif alignmethod=="min":
        alignclass = alignMin
    elif alignmethod=="centroid":
        alignclass = alignCentroid
    else:
        alignclass = alignMax
    pad = not crop
    onraw = True

    # Reference dataset
    reference = [s for s in stacks if s.endswith(refdataset)]
    if len(reference)==0:
        reference = [s for s in stacks if refdataset in s]
    if len(reference)==0:
        raise ValueError("Reference dataset \"{}\" doesn't exist.".format(refdataset))
    elif len(reference)!=1:
        # e.g. reference = Fe-K
        #      stacks = ['/detectorsum/wFe-K', '/detectorsum/Fe-K', ...]
        reference2 = [s for s in stacks if s.split('/')[-1]==refdataset]
        if len(reference2)==1:
            reference = reference2
        else:
            raise ValueError("Reference dataset \"{}\" name not specific enough (matches: {})".format(refdataset,reference))
    reference = stacks.index(reference[0])

    # Open source and destination
    bsamefile = filein==fileout
    if bsamefile:
        fin = nexus.File(filein,mode='r+')
        fout = fin
        extension = "align"

        # no need to copy groups (same file)
    else:
        fin = nexus.File(filein,mode='r')
        fout = nexus.File(fileout,mode='w' if overwrite else 'a')
        extension = ""

        if copygroups is not None:
            for grp in copygroups:
                fin.copy(fin[grp],fout)

    # Input dataset
    datasets = [fin[name][fin[name].attrs["signal"]] for name in stacks]

    # Create NXdata groups
    alignedstacks = [nexus.newNXdata(fout,s,extension) for s in stacks]

    # Output dataset names (don't exist yet)
    destlist = [grp.name+"/"+grp.attrs["signal"] for grp in alignedstacks]

    # Align
    o = alignclass(datasets,None,fout,destlist,"",stackdim=stackdim,plot=plot)
    o.align(reference,onraw = onraw,pad = pad,crop = crop,roi = roi,refimageindex=refimageindex)

    # Data sets are the default in their NXdata group
    for s in destlist:
        fout[s].attrs["signal"]=1

    # New axes
    if bsamefile:
        axesdata = [fin[a["fullname"]] for a in axes]
    else:
        axesdata = [fin[a["fullname"]][:] for a in axes]
    o.getaxesaftertransformation(axesdata)
    alignedaxes = nexus.newaxes(fout,axes,axesdata,extension)

    # Link new axes in new NXdata groups
    nexus.linkaxes(fout,alignedaxes,alignedstacks)
    alignedstack = [grp.name for grp in alignedstacks]

    # Add processing info
    if info is not None:
        info["changeofframes"] = o.absolute_cofs()
        if bsamefile:
            nexus.addinfogroup(fout,"align",info)
        else:
            nexus.copyaddinfogroup(fin,fout,"align",info)

    # Close hdf5 files
    fin.close()
    if not bsamefile:
        fout.close()

    # Return stacks and axes
    return alignedstack, alignedaxes

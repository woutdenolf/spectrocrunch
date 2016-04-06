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

from spectrocrunch.align.alignElastix import alignElastix
from spectrocrunch.align.alignSift import alignSift
from spectrocrunch.align.alignFFT import alignFFT
from . import preparenexus_hdf5_imagestacks as nexus

import h5py

def align_hdf5_imagestacks(filein,stacks,axes,fileout,alignmethod,refdataset,refimageindex=None,overwrite=False,info=None,copygroups=None):

    # Alignment method
    if alignmethod=="sift":
        alignclass = alignSift
    elif alignmethod=="elastix":
        alignclass = alignElastix
    else:
        alignclass = alignFFT
    extend = True
    onraw = True

    # Reference dataset
    reference = [s for s in stacks if s.endswith(refdataset)]
    if len(reference)==0:
        reference = [s for s in stacks if refdataset in s]

    if len(reference)==0:
        raise ValueError("Reference dataset doesn't exist.")
    elif len(reference)!=1:
        raise ValueError("Reference dataset name not specific enough.")
    reference = stacks.index(reference[0])

    # Open source and destination
    bsamefile = filein==fileout
    if bsamefile:
        fin = h5py.File(filein,'r+')
        fout = fin
        extension = "align"

        if info is not None:
            nexus.addinfogroup(fout,"align",info)
    else:
        fin = h5py.File(filein,'r')
        fout = h5py.File(fileout,'w' if overwrite else 'a')
        extension = ""

        if info is not None:
            nexus.copyaddinfogroup(fin,fout,"align",info)

        if copygroups is not None:
            for grp in copygroups:
                fin.copy(fin[grp],fout)

    # Input dataset
    datasets = [fin[name][fin[name].attrs["signal"]] for name in stacks]

    # Create NXdata groups
    alignedstacks = nexus.newgroups(fout,stacks,extension)

    # Output dataset names (don't exist yet)
    destlist = [grp.name+"/"+grp.attrs["signal"] for grp in alignedstacks]

    # Align
    o = alignclass(datasets,None,fout,destlist,"",stackdim=2)
    o.align(reference,onraw = onraw,extend = extend,refimageindex=refimageindex)

    # New axes
    if bsamefile:
        axesdata = [fin[a["fullname"]] for a in axes]
    else:
        axesdata = [fin[a["fullname"]][:] for a in axes]
    o.getaxesafteralignment(axesdata)
    alignedaxes = nexus.newaxes(fout,axes,axesdata,extension)

    # Link new axes in new NXdata groups
    nexus.linkaxes(fout,alignedaxes,alignedstacks)
    alignedstack = [grp.name for grp in alignedstacks]

    fin.close()
    if not bsamefile:
        fout.close()

    return alignedstack, alignedaxes

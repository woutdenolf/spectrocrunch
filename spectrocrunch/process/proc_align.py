# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
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

import logging

from ..h5stacks.align_hdf5_imagestacks import align_hdf5_imagestacks as alignstacks

from . import proc_common

def execute(file_in, stacks_in, axes_in, copygroups, bsamefile, default,\
            alignmethod, alignreference, refimageindex, crop, roi, plot, stackdim):

    logger = logging.getLogger(__name__)
    logger.info("Aligning image stacks ...")

    # Output file
    if bsamefile:
        file_out = file_in
    else:
        base, ext = proc_common.hdf5base(file_in)
        file_out = base+".align"+ext

    # Processing info
    info = {"method":alignmethod,\
            "pairwise":refimageindex==None,\
            "reference set":alignreference,\
            "reference image":refimageindex,\
            "crop":crop,\
            "roi":roi}

    # Align
    stacks_out, axes_out = alignstacks(file_in,stacks_in,axes_in,stackdim,file_out,alignmethod,
                                        alignreference,refimageindex=refimageindex,overwrite=True,crop=crop,
                                        roi=roi,plot=plot,info=info,copygroups=copygroups)

    # Default
    proc_common.defaultstack(file_out,stacks_out,default)
    
    return file_out, stacks_out, axes_out


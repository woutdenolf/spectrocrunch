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

from spectrocrunch.h5stacks.math_hdf5_imagestacks import crop_hdf5_imagestacks as cropstacks

from . import proc_common

def execute(file_in, stacks_in, axes_in, copygroups, bsamefile, default,\
            cropinfo):

    logger = logging.getLogger(__name__)
    logger.info("Cropping image stacks ...")

    # Output file
    if bsamefile:
        file_out = file_in
    else:
        base, ext = proc_common.hdf5base(file_in)
        file_out = base+".crop"+ext

    # Processing info
    if "nanval" in cropinfo:
        info = {"crop value":cropinfo["nanval"],"reference set":cropinfo["reference set"]}
    else:
        info = {"image roi":cropinfo["roi"]}

    # Align
    stacks_out, axes_out = cropstacks(file_in,file_out,axes_in,stacks_in,stacks_in,cropinfo,overwrite=True,info=info,copygroups=copygroups)

    # Default
    proc_common.defaultstack(file_out,stacks_out,default)
    
    return file_out, stacks_out, axes_out


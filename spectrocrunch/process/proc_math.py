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

from ..h5stacks.math_hdf5_imagestacks import math_hdf5_imagestacks as math

from . import proc_common

def execute(file_in,stacks_in,axes_in,copygroups,bsamefile,default,\
            expression,skipnames,stackdim=None,copyskipped=True,extension="math"):
    logger = logging.getLogger(__name__)
    logger.info("Normalization ...")

    # Stacks to normalize
    innames = proc_common.selectnotgroups(stacks_in,skipnames)
    if len(innames)==0:
        raise ValueError("No other stacks found than these: "+str.join(", ",skipnames))

    # Output file
    if bsamefile:
        file_out = file_in
    else:
        base, ext = proc_common.hdf5base(file_in)
        file_out = base+"."+extension+ext

    # Expand full paths
    fixednames = math.extractexpression(expression,allowempty=False)
    fixedpaths = proc_common.selectgroups(stacks_in,fixednames,sort=True)
    expression = math.replaceexpression(expression,fixednames,fixedpaths)

    # Processing info
    info = {}
    info["expression"] = expression

    # Normalize
    stacks_out, axes_out = math.calc_hdf5_imagestacks(file_in,file_out,axes_in,expression,innames,innames,overwrite=True,info=info,copygroups=copygroups,stackdim=stackdim,extension=extension)

    # Copy unnormalized stacks when new file
    if copyskipped:
        innames = proc_common.selectgroups(stacks_in,skipnames)
        if len(innames)!=0:
            if file_in==file_out:
                stacks_out += innames
            else:
                tmp_stacks, tmp = math.copy_hdf5_imagestacks(file_in,file_out,axes_in,innames,innames,overwrite=False)
                stacks_out += tmp_stacks

    # Default
    proc_common.defaultstack(file_out,stacks_out,default)
    
    return file_out, stacks_out, axes_out


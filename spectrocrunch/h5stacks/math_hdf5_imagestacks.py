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

import math_hdf5_imagestacks_expression as m_expression
import math_hdf5_imagestacks_copy as m_copy
import math_hdf5_imagestacks_crop as m_crop
import math_hdf5_imagestacks_replace as m_replace
import math_hdf5_imagestacks_resample as m_resample

from ..io import nexus as nexus
from ..common.integerbase import integerbase

from spectrocrunch.common.Enum import Enum
operationType = Enum(['expression','copy','crop','replace','resample'])

import re
import logging

def math_hdf5_imagestacks(filein,fileout,axes,operation,varargs,fixedargs,ret,extension="",overwrite=False,info=None,copygroups=None):
    """Perform some operations on hdf5 imagestacks

    Args:
        filein:
        fileout:
        axes(list(str)): list of NXdata groups
        operation(dict): type, value
        varargs(list(dict)|list(array-like)): 
        fixedargs(dict):
        ret: paths names of the results (same elements as varargs)
        extension:
        overwrite:
        info:
        copygroups:
    Returns:
        tuple(list,list)
    """

    bsame = filein==fileout

    # Open files
    if bsame:
        fin = nexus.File(filein,mode='r+')
        fout = fin
    
        if info is not None:
            nexus.addinfogroup(fout,"math",info)
    else:
        fin = nexus.File(filein,mode='r')
        fout = nexus.File(fileout,mode='w' if overwrite else 'a')
        extension = ""

        if info is not None:
            nexus.copyaddinfogroup(fin,fout,"math",info)

        if copygroups is not None:
            for grp in copygroups:
                fin.copy(fin[grp],fout)

    # Create NXdata groups
    retstacks = [nexus.newNXdata(fout,s,extension) for s in ret]

    # Operation
    axesdata = []
    if operation["type"]==operationType.expression:
        m_expression.evaluate(operation,fin,varargs,fixedargs,retstacks)
    elif operation["type"]==operationType.crop:
        axesdata = m_crop.evaluate(operation,fin,varargs,retstacks,axes)
    elif operation["type"]==operationType.copy:
        m_copy.evaluate(operation,fin,varargs,retstacks)
    elif operation["type"]==operationType.replace:
        m_replace.evaluate(operation,fin,varargs,retstacks)
    elif operation["type"]==operationType.resample:
        axesdata = m_resample.evaluate(operation,fin,varargs,retstacks,axes)

    # New axes
    if len(axesdata)==0:
        # Copy of the old axes
        if bsame:
            axesdata = [fin[a["fullname"]] for a in axes]
        else:
            axesdata = [fin[a["fullname"]][:] for a in axes]

    retaxes = nexus.newaxes(fout,axes,axesdata,extension)

    # Link groups to axes
    nexus.linkaxes(fout,retaxes,retstacks)

    # Keep only names
    retstacks = [dset.name for dset in retstacks]

    # Return
    fin.close()
    if not bsame:
        fout.close()

    return retstacks, retaxes

def extractexpression(expression,allowempty=True):
    if allowempty:
        return list(set(re.findall("\{(.*?)\}",expression)))
    else:
        return list(set(re.findall("\{([^\{\}]+)\}",expression)))

def replaceexpression(expression,olds,news):
    for old,new in zip(olds,news):
        expression = expression.replace("{{{}}}".format(old),"{{{}}}".format(new))
    return expression

def parseexpression(expression,stacks):
    """
    Args:
        expression(str): the variable argument is indicated by "{}" 
        fixed arguments are indeicate by "{name}"
        e.g. "{}*nanone({icr}/({I0}*{ocr})"
        stacks(list(str)): list of variable arguments, e.g. ["icr","ocr","I0",...]
    """
    _vars = extractexpression(expression)
    if "" not in _vars:
        raise ValueError("The expression does not contain {} which indicates the variable argument.")

    digs = ['b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']
    o = integerbase(digs = digs)
    
    # Variable argument
    expression = expression.replace("{}","var_a") # do not use format
    _vars.remove("")
    n = len(stacks)
    varargs = [None]*n
    for i in range(n):
        varargs[i] = {"a":stacks[i]}

    # Fixed arguments
    fixedargs = {}
    for i in range(len(_vars)):
        varname = o.int2base(i+1)
        expression = expression.replace("{{{}}}".format(_vars[i]),"var_{}".format(varname))
        fixedargs[varname] = _vars[i]

    return expression,varargs,fixedargs

def calc_hdf5_imagestacks(filein,fileout,axes,expression,stacks,retstack,overwrite=False,info=None,copygroups=None,stackdim=None,extension=""):
    logger = logging.getLogger(__name__)
    logger.debug("Stack math expression: {}".format(expression))

    expression,varargs,fixedargs = parseexpression(expression,stacks)
    operation = {"type":operationType.expression,"value":expression,"stackdim":stackdim,"sliced":stackdim is not None}
    return math_hdf5_imagestacks(filein,fileout,axes,operation,varargs,fixedargs,retstack,extension=extension,overwrite=overwrite,info=info,copygroups=copygroups)

def minlog_hdf5_imagestacks(filein,fileout,axes,stacks,retstack,overwrite=False,info=None,copygroups=None,stackdim=None):
    expression = "-ln(var_a)"
    operation = {"type":operationType.expression,"value":expression,"stackdim":stackdim,"sliced":stackdim is not None}

    nI = len(stacks)
    varargs = [None]*nI
    for i in range(nI):
        varargs[i] = {"a":stacks[i]}

    return math_hdf5_imagestacks(filein,fileout,axes,operation,varargs,{},retstack,extension="minlog",overwrite=overwrite,info=info,copygroups=copygroups)

def copy_hdf5_imagestacks(filein,fileout,axes,stacks,retstack,overwrite=False,info=None,copygroups=None,stackdim=None):
    operation = {"type":operationType.copy,"sliced":stackdim is not None,"stackdim":stackdim}
    return math_hdf5_imagestacks(filein,fileout,axes,operation,stacks,{},retstack,extension="copy",overwrite=overwrite,info=info,copygroups=copygroups)

def crop_hdf5_imagestacks(filein,fileout,axes,stacks,retstack,cropinfo,overwrite=False,info=None,copygroups=None):
    operation = cropinfo
    operation["type"] = operationType.crop
    return math_hdf5_imagestacks(filein,fileout,axes,operation,stacks,{},retstack,extension="crop",overwrite=overwrite,info=info,copygroups=copygroups)

def replacevalue_hdf5_imagestacks(filein,fileout,axes,stacks,retstack,orgvalue,newvalue,overwrite=False,info=None,copygroups=None,stackdim=None):
    operation = {"type":operationType.replace,"v1":orgvalue,"v2":newvalue,"sliced":stackdim is not None,"stackdim":stackdim}
    return math_hdf5_imagestacks(filein,fileout,axes,operation,stacks,{},retstack,extension="replace",overwrite=overwrite,info=info,copygroups=copygroups)

def resample_hdf5_imagestacks(filein,fileout,axes,stacks,retstack,resampleinfo,overwrite=False,info=None,copygroups=None):
    operation = resampleinfo
    operation["type"] = operationType.resample
    return math_hdf5_imagestacks(filein,fileout,axes,operation,stacks,{},retstack,extension="resample",overwrite=overwrite,info=info,copygroups=copygroups)




# -*- coding: utf-8 -*-
#
#   Copyright (C) 2017 European Synchrotron Radiation Facility, Grenoble, France
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

import random
import numpy as np
import numbers
import itertools

def isadvanced(index):
    """Check for advanced indexing

    Args:
        index(index): object used in indexing or slicing
    Returns:
        num: length of the list indexing
    """

    if isinstance(index,tuple):
        return max([isadvanced(i) for i in index])
    elif isinstance(index,list):
        return len(index)
    else:
        return 0

def nadvanced(index):
    """Number of dimensions with advanced indexing

    Args:
        index(index): object used in indexing or slicing
    Returns:
        num: length of the list indexing
    """

    if isinstance(index,tuple):
        return sum([nadvanced(i) for i in index])
    elif isinstance(index,list):
        return 1
    else:
        return 0

def adjacent(arr):
    """Check whether list indexing is done one adjacent dimensions or not

        >>> a=np.zeros((10,20,30,40))
        >>> a[[0,1],:,:,[0,1]].shape == (2, 20, 30) # not adjacent
        >>> a[:,[0,1],[0,1],:].shape == (10, 2, 40) # adjacent

    Args:
        arr(array): index or bool array
    Returns:
        bool
    """
    if isinstance(arr,tuple):
        arr = [isinstance(i,list) for i in arr]

    return sum([k for k,g in itertools.groupby(arr)])==1

def nonchanging(index,shape=None):
    """Check whether slicing takes the entire range or a subrange

    Args:
        index(index): object used in indexing or slicing
        shape(Optional(tuple)): dimension before applying index
    Returns:
        bool
    """
    if isinstance(index,tuple):
        if nadvanced(index)>1:
            return False

        if shape is None:
            return all(nonchanging(i) for i in index)
        else:
            return all(nonchanging(i,s) for i,s in zip(index,shape))
    elif isinstance(index,slice):
        return index==slice(None) or index==slice(0,None) or index==slice(None,None,1)
    elif index is Ellipsis: # slice(None) for missing dimensions in tuple index
        return True
    elif index is np.newaxis: # adds a dimension
        return False
    elif isinstance(index,list):
        if shape is None:
            return False # could be True, but we can't know when shape is not given
        else:
            return index.count(True)==len(index) or index==range(shape)
    else:
        return False

def positiveaxis(axis,ndim):
    """Positive axis

    Args:
        axis(num): dimension index
        ndim(num): number of dimensions
    Returns:
        num
    """
    if axis<0:
        axis += ndim
    if axis<0 or axis>=ndim:
        raise IndexError("axis out of range")
    return axis
    
def expand(index,ndim):
    """Expand index to ndim dimensions

    Args:
        index(index): object used in indexing or slicing
        ndim(num): number of dimensions
    Returns:
        index
    """
    if isinstance(index,tuple):
        nindex = sum([i is not np.newaxis for i in index])
        if nindex>ndim:
            raise IndexError("invalid index")

        if Ellipsis in index:
            i = index.index(Ellipsis)
            index = index[:i] + (slice(None),)*(ndim-nindex+1) + index[i+1:]
        elif nindex!=ndim:
            index = index + (slice(None),)*(ndim-nindex)

    elif index is Ellipsis:
        index = (slice(None),)*ndim

    elif index is np.newaxis:
        index = (None,)+(slice(None),)*ndim

    elif ndim>1:
        index = (index,)+(slice(None),)*(ndim-1)
    
    return index

def axisindex(index,axis,ndim):
    """Find axis in expanded index

    Args:
        axis(num): dimension index
        ndim(num): number of dimensions
    Returns:
        num
    """

    axis = positiveaxis(axis,ndim)
    lst = list(np.cumsum([i is not np.newaxis for i in index]))
    return lst.index(axis+1)

def fulldim(index,ndim):
    """Dimension after indexing (without squeezing singletons)

    Args:
        index(index):
        ndim(num): number of dimensions
    Returns:
        num
    """
    return ndim + sum([i is None for i in index])

def replace(index,ndim,axis,rindex):
    """Replace indexing for a specified dimension

    Args:
        index(index): object used in indexing or slicing
        ndim(num): number of dimensions
        axis(num): dimension to be replaced
        rindex(slice, num, ...): new indexing for this dimensions
    Returns:
        index
    """

    index2 = expand(index,ndim)

    axis = axisindex(index2,axis,ndim)

    return index2[:axis] + (rindex,) + index2[axis+1:]

def replacefull(index,ndim,axis):
    """Set the specified dimension to full range

    Args:
        index(index): object used in indexing or slicing
        ndim(num): number of dimensions
        axis(num): dimension to be set to full range
    Returns:
        index
    """

    return replace(index,ndim,axis,slice(None))

def extract_advanced(index):
    """Replace advanced indexing

        >>> index = (slice(1,2),range(5),slice(3,4),range(5))
        >>> index1 = tuple([i if isinstance(i,list) else slice(None) for i in index])
        >>> index2,_ = extract_advanced(index)
        >>> a[index] == a[index1][index2]

    Args:
        index(index): object used in indexing or slicing (expanded) 
    Returns:
        index
    """

    i = None

    if isinstance(index,tuple):
        b = [isinstance(j,list) for j in index]
        nlists = sum(b)
        if nlists>0:
            if adjacent(b):
                i = b.index(True)
            else:
                i = 0
            index = [ind for j,ind in enumerate(index) if not b[j]]
            index.insert(i,slice(None))
            index = tuple(index)

    return index,i

def advanced_index_axis(index):
    """Get list-axis after applying advanced indexing

    Args:
        index(index): object used in indexing or slicing (expanded) 
    Returns:
        num
    """
    i = None

    if isinstance(index,tuple):
        b = [isinstance(j,list) for j in index]
        if any(b):
            if adjacent(b):
                i = b.index(True)
            else:
                i = 0
    elif isinstance(index,list):
        i = 0

    return i

def unsqueeze_index(index,ndim):
    """Restore singleton dimensions after they have been squeezed.
       
    Args:
        index(index): object used in indexing or slicing
        ndim(num): number of dimensions
    Returns:
        index or None
    """

    # Expanded index
    index = expand(index,ndim)

    # More than one list
    index,_ = extract_advanced(index)

    # Add new axis for singleton indexing
    index = [np.newaxis if isinstance(i,numbers.Number) else slice(None) for i in index]

    return tuple(index)

def unsqueeze(arr,index,ndim):
    """Restore singleton dimensions after they have been squeezed.

    Args:
        arr(array): array after applying index
        index(index): object used in indexing or slicing
        ndim(num): number of dimensions before indexing
    Returns:
        array
    """
    return arr[unsqueeze_index(index,ndim)]

def expanddims(arr,ndim):
    """Expand array dimensions

    Args:
        arr(array):
        ndim(num): new dimensions
    Returns:
        array
    """
    n = arr.ndim-ndim
    if n > 1:
        arr = arr[(Ellipsis,)+(np.newaxis,)*n]
    return arr

class op_index(object):
    def __init__(self,index=None):
        self.index = index

    def __call__(self,data):
        if self.index is not None:
            data = data[self.index]
        return data

    def __str__(self):
        return "indexing: {}".format(self.index)

class op_transpose(object):

    def __init__(self,axes=None):
        self.axes = axes

    def __call__(self,data):
        if self.axes is not None:
            data = np.transpose(data,self.axes)

        return data

    def __str__(self):
        return "transpose: {}".format(self.axes)

class operators(object):

    def __init__(self):
        self.ops = []

    def append(self,op):
        if isinstance(op,list):
            self.ops += op
        elif isinstance(op,operators):
            self.ops += op.ops
        else:
            self.ops.append(op)

    def __call__(self,data):
        for o in self.ops:
            data = o(data)

        return data

def move_axis(ndim,i1,i2):
    """ Transpose axes to move axis i1 to i2

    Args:
        ndim(num):
        i1(num):
        i2(num):

    Returns:
        list
    """
    if i1==i2:
        return None
    axes = range(ndim)

    axes.insert(i2,axes.pop(i1))
    return axes
    
def extract_newaxis(index):
    """Extract np.newaxis from index and applied later.

        >>> index2,addnewaxis = extract_newaxis(index1)
        >>> unsqueeze(a[index1],index1,a.ndim) == addnewaxis(unsqueeze(a[index2],index2,a.ndim))

    Args:
        index(index): object used in indexing or slicing (expanded)
    Returns:
        2-tuple
    """

    ops = operators()

    bnewaxes = any(i is np.newaxis for i in index)
    if bnewaxes:
        # index after applying the advanced part of the index
        indexadv,iadv = extract_advanced(index)

        # (1) index without np.newaxis
        index = tuple([i for i in index if i is not np.newaxis])
        iadv0 = advanced_index_axis(index)

        # (2) add new dimensions
        index_newaxes = tuple([i if i is np.newaxis else slice(None) for i in indexadv])
        ops.append(op_index(index_newaxes))
        if iadv0 is not None:
            iadv0 += sum([i is np.newaxis for i in index_newaxes[:iadv0+1]])

        # (3) shift destination axis for advanced indexing
        axes = move_axis(len(indexadv),iadv0,iadv)
        if axes is not None:
            ops.append(op_transpose(axes))

    return index,ops

def popaxis(shape,axis):
    shape = list(shape)
    shapeaxis = shape.pop(axis)
    shape = tuple(shape)
    return shapeaxis,shape

def extract_axis(index,ndim,axis,shapefull):
    """Extract axis from index and applied in stacking later.

    Args:
        index(index): object used in indexing or slicing (expanded and no np.newaxis)
        ndim(num): number of dimensions before indexing
        axis(num): axis index (positive)
        shapefull(tuple): full stack dimension (without applying index)
    Returns:
        tuple
    """

    ops = operators()

    # index of the list-axis after advanced indexing
    indexadv,iadv = extract_advanced(index)

    # separate axis from other dimensions
    indexaxis,indexrest = popaxis(index,axis)
    if shapefull is None:
        shapeaxis = None
        shaperest = None
    else:
        shapeaxis,shaperest = popaxis(shapefull,axis)

    if isadvanced(indexaxis):
        # axis dimension is advanced and will the advanced destination dimension
        axis = iadv
    else:
        # axis dimension is not advanced
        if iadv is not None:
            # indexrest is advanced

            # Extracting axis can cause non-adjacent to be adjacent
            if not adjacent(index) and adjacent(indexrest):
                axes = move_axis(len(indexadv),advanced_index_axis(indexrest),iadv)
                if axes is not None:
                    ops.append(op_transpose(axes))

            # Number of advanced indices before axis
            ndisappear = sum([isinstance(i,list) for i in index[:axis]])

            # Advanced destination dimension before axis
            ndisappear -= iadv<=axis

            axis -= ndisappear

    return indexrest,shaperest,indexaxis,shapeaxis,axis,ops

def prepare_slicedstack(index,ndim,axis,shapefull):
    """

    Args:
        index(index): object used in indexing or slicing
        ndim(num): number of dimensions before indexing
        axis(num): generator axis
        shapefull(tuple): full stack dimension (without applying index)
    Returns:
        tuple
    """

    index = expand(index,ndim)

    axis = positiveaxis(axis,ndim)

    index,ops2 = extract_newaxis(index)

    indexrest,shaperest,indexaxis,shapeaxis,axis,ops1 = extract_axis(index,ndim,axis,shapefull)

    ops1.append(ops2)

    return indexrest,shaperest,indexaxis,shapeaxis,axis,ops1

def slicedstack(generator,arglist,index,ndim,shapefull=None,axis=-1):
    """Similar to np.stack(...,axis=...)[index] except that the data
       still needs to be generate by looping over the axis dimension.
       The data generation is minimized based in index.

       >>> data = np.stack(map(generator,args),axis=-1)
       >>> shapefull = data.shape
       >>> ndim = data.ndim
       >>> data2 = slicedstack(generator,args,index,ndim,shapefull=shapefull,axis=-1)
       >>> np.testing.assert_array_equal(indexing.unsqueeze(data[index],index,ndim),data2)

    Args:
        generator: data generator
        arglist(list): generator arguments
        index(index): object used in indexing or slicing
        ndim(num): stack dimensions (without applying index)
        shapefull(Optional(tuple)): full stack dimension (without applying index)
        axis(Optional(num)): axis corresponding to arglist

    Returns:
        array
    """
    
    # Split axis indexing from indexing other dimensions
    indexrest,shaperest,indexaxis,shapeaxis,axis,ops = prepare_slicedstack(index,ndim,axis,shapefull)
    ndimrest = ndim-1

    # Apply indexaxis on axis dimension
    if not nonchanging(indexaxis,shape=shapeaxis):
        if isinstance(indexaxis,list):
            # No advanced indexing for list so use comprehension
            arglist = [arglist[i] for i in indexaxis]
        elif isinstance(indexaxis,numbers.Number):
            # No singleton dimension squeezing
            arglist = [arglist[indexaxis]]
        else:
            arglist = arglist[indexaxis]

    # Apply indexrest on the other dimensions 
    if nonchanging(indexrest,shape=shaperest):
        data = [expanddims(generator(arg),ndimrest) for arg in arglist]
    else:
        # Handle advanced indexing
        na = isadvanced(indexaxis)
        if na!=0:
            squeezeadvanced = tuple([0 if isinstance(i,list) else slice(None) for i in indexrest])

            indexlist = [tuple([i[j] if isinstance(i,list) else i for i in indexrest]) for j in range(na)]

            data = [unsqueeze(expanddims(generator(arg),ndimrest)[ind],ind,ndimrest)[squeezeadvanced] for arg,ind in zip(arglist,indexlist)]
        else:
            data = [unsqueeze(expanddims(generator(arg),ndimrest)[indexrest],indexrest,ndimrest) for arg in arglist]

    # Stack
    if len(data)==0:
        tmp,_ = extract_advanced(expand(index,ndim))

        ndim = len(tmp)
        shape = [0]*ndim
        data = np.empty(tuple(shape))
    else:
        data = ops(np.stack(data,axis=axis))

    return data


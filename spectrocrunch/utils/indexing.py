import numpy as np
import numbers

from . import instance
from . import listtools


def isadvanced(index):
    """Check for advanced indexing

    Args:
        index(index): object used in indexing or slicing

    Returns:
        bool
    """
    if isinstance(index, tuple):
        return any(isadvanced(ind) for ind in index)
    return instance.islistgen(index)


def lengthadvanced(index, shape):
    """Length after advanced indexing

    Args:
        index(index): object used in indexing or slicing (expanded and no np.newaxis)
        shape: dimensions

    Returns:
        num:
    """
    if isinstance(index, tuple):
        tmp = [ind for ind in index if ind is not np.newaxis]
        tmp = [lengthadvanced(ind, s) for ind, s in zip(tmp, shape)]
        if not any(tmp):
            return 0
        tmp = [t for t in tmp if t != 0]
        if tmp.count(tmp[0]) != len(tmp):
            raise IndexError(
                "shape mismatch: indexing arrays could not be broadcast together with shapes {}".format(
                    tmp
                )
            )
        return tmp[0]
    elif instance.islistgen(index):
        if instance.isboolsequence(index):
            return sum(index)
        else:
            return len(index)
    else:
        return 0


def indexcount(index, tp):
    """Number of dimensions with a specific index type

    Args:
        index(index): object used in indexing or slicing
    Returns:
        num: length of the list indexing
    """
    if isinstance(index, tuple):
        return sum([indexcount(i, tp) for i in index])
    elif isinstance(index, tp):
        return 1
    else:
        return 0


def nadvanced(index):
    """Number of dimensions with advanced indexing

    Args:
        index(index): object used in indexing or slicing
    Returns:
        num: length of the list indexing
    """
    return indexcount(index, instance.listgentypes)


def nsingleton(index):
    """Number of singleton dimensions

    Args:
        index(index): object used in indexing or slicing
    Returns:
        num: length of the list indexing
    """
    return indexcount(index, numbers.Number)


def isdimchanging(index):
    """Check whether a dimension may be squeezed after indexing

    Args:
        index(index): object used in indexing or slicing
    Returns:
        array or num
    """
    if isinstance(index, tuple):
        return [isdimchanging(ind) for ind in index]
    else:
        return (
            instance.islistgen(index)
            or isinstance(index, numbers.Number)
            or index is np.newaxis
        )


def replace_dimnonchanging(index):
    if isinstance(index, tuple):
        return tuple([replace_dimnonchanging(ind) for ind in index])
    elif isdimchanging(index):
        return index
    else:
        return slice(None)


def extract_dimchanging(index):
    """Separate indices that preserve the dimensions from indexing that don't

        >>> index = (slice(1,2),range(5),slice(3,4),range(5))
        >>> index1 = extract_dimchanging(index)
        >>> index2,iadv = extract_dimnonchanging(index)
        >>> a[index] == a[index1][index2]

    Args:
        index(index): object used in indexing or slicing (expanded)
    Returns:
        index: index1
    """
    return replace_dimnonchanging(index)


def replace_dimchanging(index):
    if isinstance(index, tuple):
        return tuple([replace_dimchanging(ind) for ind in index])
    elif isadvanced(index):
        return [0, 0], slice(None), False
    elif isinstance(index, numbers.Number):
        return 0, slice(None), False
    elif index is np.newaxis:
        return np.newaxis, slice(None), True
    else:
        return slice(None), index, True


def extract_dimnonchanging(index):
    """Separate indices that preserve the dimensions from indexing that don't

        >>> index = (slice(1,2),range(5),slice(3,4),range(5))
        >>> index1 = extract_dimchanging(index)
        >>> index2,iadv = extract_dimnonchanging(index)
        >>> a[index] == a[index1][index2]

    Args:
        index(index): object used in indexing or slicing (expanded)
    Returns:
        index,num: index2, final axis of list indexing
    """
    if len(index) == 0:
        return (), None

    index1, index2, bkeepindex2 = zip(*replace_dimchanging(index))
    index2 = listtools.listadvanced_bool(index2, bkeepindex2)

    s = np.empty((1,) * len(index1))[index1].shape

    # 2 -> destination axis of advanced indexing
    # 1 -> any other dimension which is not squeezed
    if 2 in s:
        iadv = s.index(2)
    else:
        iadv = None

    if iadv is not None:
        index2.insert(iadv, slice(None))

    return tuple(index2), iadv


def replace_nonnewaxis(index):
    if isinstance(index, tuple):
        return tuple([replace_nonnewaxis(ind) for ind in index])
    elif isadvanced(index):
        return [0, 0], [0, 0], True
    elif isinstance(index, numbers.Number):
        return 0, 0, True
    elif index is np.newaxis:
        return np.newaxis, np.newaxis, False
    else:
        return slice(None), slice(None), True


def extract_newaxis(index):
    """Extract np.newaxis from index and applied later.

        >>> index2,addnewaxis = extract_newaxis(index1)
        >>> a[index1] == addnewaxis(a[index2])

    Args:
        index(index): object used in indexing or slicing (expanded)
    Returns:
        2-tuple
    """
    ops = operators()
    bnew = [ind is np.newaxis for ind in index]
    if any(bnew):
        indexwithnew, indexwithoutnew, bindexwithoutnew = zip(
            *replace_nonnewaxis(index)
        )
        indexwithoutnew = tuple(
            listtools.listadvanced_bool(indexwithoutnew, bindexwithoutnew)
        )

        s = (3,) * len(indexwithoutnew)
        s1 = np.empty(s)[indexwithnew].shape
        # 1 -> new dimension
        # 2 -> destination axis of advanced indexing
        # 3 -> any other dimension which is not squeezed

        s1b = [i for i in s1 if i != 1]
        if len(s) == 0:
            s2 = ()
        else:
            s2 = np.empty(s)[indexwithoutnew].shape
        # 2 -> destination axis of advanced indexing
        # 3 -> any other dimension which is not squeezed

        # (1) Index without np.newaxis
        index = tuple(listtools.listadvanced_bool(index, bnew, bnot=True))

        # (2) Transpose
        if 2 in s2:
            axes = move_axis(len(s2), s2.index(2), s1b.index(2))
            if axes is not None:
                ops.append(op_transpose(axes))

        # (3) Add new dimensions
        index_newaxes = tuple([np.newaxis if i == 1 else slice(None) for i in s1])
        ops.append(op_index(index_newaxes))
    return index, ops


def popaxis(shape, axis):
    shape = list(shape)
    shapeaxis = shape.pop(axis)
    shape = tuple(shape)
    return shapeaxis, shape


def replace_axesorder(index):
    if isinstance(index, tuple):
        return tuple([replace_axesorder(ind) for ind in index])
    elif isadvanced(index):
        return []
    elif isinstance(index, numbers.Number):
        return 0
    elif index is np.newaxis:
        return np.newaxis
    else:
        return slice(None)


def axesorder_afterindexing(index, ndim):
    """Axes order after indexing specified (negative axes are np.newaxis).

    Args:
        index(index): object used in indexing or slicing (expanded and no np.newaxis)
        ndim(num): dimensions before indexing

    Returns:
        list:
        index:
    """
    indexexpanded = expand(index, ndim)

    # Dimensions which are indexed by lists
    indexnonew = tuple([ind for ind in indexexpanded if ind is not np.newaxis])
    ilist = [i for i, ind in enumerate(indexnonew) if isadvanced(ind)]

    # Dimensions are recognized by their length:
    # 0     -> destination of advanced indexing
    # 1     -> np.newaxis
    # else  -> original dimension (2 is the first)
    shape = tuple(range(2, ndim + 2))
    arr = np.empty(shape)
    shape = arr[replace_axesorder(indexexpanded)].shape

    # Axis origins
    axes = list(shape)
    j = -1
    for i, n in enumerate(axes):
        if n == 0:
            axes[i] = ilist
        elif n == 1:
            axes[i] = j
            j -= 1
        else:
            axes[i] = n - 2

    return axes, indexnonew


def shape_afterindexing(shape, index, ndim=None):
    """Shape after indexing

    Args:
        shape(tuple):
        index(index):
        ndim(Optional(int)): in case shape is None

    Returns:
        tuple
    """
    if shape is None:
        shape = (1,) * ndim

    axes, indexnonew = axesorder_afterindexing(index, len(shape))

    # Expected shape after indexing
    s2 = [0] * len(axes)
    for i, iaxes in enumerate(axes):
        if instance.isnumber(iaxes) and iaxes < 0:
            s2[i] = 1
        elif instance.islistgen(iaxes):
            s2[i] = lengthadvanced(indexnonew, shape)
        else:
            ind = indexnonew[iaxes]
            if isinstance(ind, slice):
                s2[i] = len(range(shape[iaxes])[ind])
            else:
                s2[i] = shape[iaxes]

    return tuple(s2)


def replace_trimindex(index):
    if isinstance(index, tuple):
        return tuple([replace_trimindex(ind) for ind in index])
    elif isadvanced(index):
        return [0, 0]
    elif isinstance(index, numbers.Number):
        return 0
    else:
        return slice(None)


def extract_axis(index, axis, shapefull):
    """Extract axis from index and applied in stacking later.

    Args:
        index(index): object used in slicing (expanded and no np.newaxis)
        axis(num): stack axis before indexing
        shapefull(tuple): full stack dimension (without applying index)

    Returns:
        indexrest(index): index dimensions != axis
        shaperest(tuple): shape dimensions != axis
        indexaxis(index): index dimensions == axis
        shapeaxis(tuple): shape dimensions == axis
        axis(num|None): stack axis after indexing
        ops(operators): transpose or identity operator
    """
    ops = operators()
    ndim = len(index)
    axis = positiveaxis(axis, ndim)

    # separate axis from other dimensions
    indexaxis, indexrest = popaxis(index, axis)
    if shapefull is None:
        shapeaxis = None
        shaperest = None
    else:
        shapeaxis, shaperest = popaxis(shapefull, axis)

    # index after dimension modifying index
    indexadv, iadv = extract_dimnonchanging(index)
    if isadvanced(indexaxis):
        # axis dimension is advanced and will be the advanced destination dimension
        axis = iadv
    else:
        # axis dimension is not advanced
        s = [1] * ndim
        s[axis] = 3
        s1 = np.empty(tuple(s))[replace_trimindex(index)].shape
        # 1 -> any other dimension
        # 2 -> destination axis of advanced indexing
        # 3 -> axis after indexing (not always there)

        if 3 in s1:  # axis not singleton
            axis = s1.index(3)
        else:
            axis = None

        if 2 in s1:
            s = [1] * (ndim - 1)
            s2 = np.empty(tuple(s))[replace_trimindex(indexrest)].shape
            # 1 -> any other dimension
            # 2 -> destination axis of advanced indexing

            if 2 in s2:
                # nstack([ .... ],axis=axis)
                if axis is not None:
                    s2 = list(s2)
                    s2.insert(axis, 3)

                # Transpose
                axes = move_axis(len(s1), s2.index(2), s1.index(2))
                if axes is not None:
                    ops.append(op_transpose(axes))

    return indexrest, shaperest, indexaxis, shapeaxis, axis, ops


class op_singletonindex(object):
    """Apply singelton index to certain axes and restore axes when requested"""

    def __init__(self, selaxes, restore):
        self.selaxes = selaxes
        self.restoreaxes = np.asarray(selaxes)[restore]

    def __str__(self):
        return "singletonindex: select axes = {}, restore axes = {}".format(
            self.selaxes, self.restoreaxes
        )

    def __call__(self, data, selind):
        n = len(self.restoreaxes)
        breshape = n != 0

        if breshape:
            shape = np.asarray(data.shape)
            b = np.full(len(shape), True)

            shape[self.selaxes] = 0  # squeeze
            b[self.selaxes] = False

            shape[self.restoreaxes] = 1  # restore
            b[self.restoreaxes] = True

            shape = shape[b]
            shape = tuple(shape)

        indexextract = np.full(data.ndim, slice(None))
        ret = []
        for ind in selind:
            if len(ind) != len(self.selaxes):
                raise RuntimeError(
                    "We have {} selected axes but {} indices.".format(
                        len(self.selaxes), len(ind)
                    )
                )

            # Index the selaxes
            indexextract[self.selaxes] = ind
            d = data[tuple(indexextract)]

            # Restore singelton dimensions that weren't already singletons
            if breshape:
                if d.size == 0:
                    d = np.empty(shape, dtype=d.dtype)
                else:
                    d = d.reshape(shape)
            ret.append(d)

        return ret


def replacefull_transform(index, fullaxes, ndim, restoreadvanced=True):
    """

    Args:
        index(index): object used in slicing
        fullaxes(list): axes to not index
        ndim(num): dimensions before indexing
    Returns:
        list,operators
    """
    # Index after replacing
    fullaxes = [positiveaxis(axis, ndim) for axis in fullaxes]
    indexfull = replacefull(index, ndim, fullaxes)

    # Dimensions with and without replacing
    axes1, _ = axesorder_afterindexing(index, ndim)
    axes2, _ = axesorder_afterindexing(indexfull, ndim)

    # Advanced indexing dimensions
    i1list = listtools.where(axes1, lambda x: instance.islistgen(x))
    if len(i1list) == 0:
        laxes1 = []
    else:
        laxes1 = axes1[i1list[0]]

    # Check for squeezed and advanced dimensions
    faxes1 = list(listtools.flatten(axes1))
    restore = [False] * len(fullaxes)
    for i, fullaxis in enumerate(fullaxes):
        if fullaxis in faxes1:
            # full dimension would have been preserved
            restore[i] = True

            if fullaxis in laxes1:
                # full dimension would have been advanced

                if laxes1 == [fullaxis]:
                    # Only advanced index: keep its place
                    restore[i] = restoreadvanced
                    axes1[i1list[0]] = fullaxis

                    laxes1 = []
                else:
                    # Not the only advanced index: extract and append
                    restore[i] = False
                    axes1[i1list[0]] = [a for a in laxes1 if a != fullaxis]
                    axes1.append(fullaxis)

                    laxes1 = axes1[i1list[0]]
        else:
            # Replaced dimension would have been squeezed
            restore[i] = False

            # Put squeezed dimension at the end
            axes1.append(fullaxis)

    # Transpose data after indexfull to match data after index
    ind = [axes2.index(a) for a in axes1 if a in axes2]
    if ind != list(range(len(ind))) and sorted(ind) == list(range(len(ind))):
        postindexfull = op_transpose(ind)
        axes2 = listtools.listadvanced_int(axes2, ind)
    else:

        def postindexfull(x):
            return x

    # For extracting particular indices along the full dimensions
    selaxes = [axes2.index(a) for a in fullaxes]
    singletonindex = op_singletonindex(selaxes, restore)
    return indexfull, postindexfull, singletonindex


def nonchanging(index, shape=None):
    """Check whether slicing is changing any axes (size or order)

    Args:
        index(index): object used in slicing
        shape(Optional(tuple)): dimension before applying index
    Returns:
        bool
    """
    if isinstance(index, tuple):
        na = nadvanced(index)
        ns = nsingleton(index)
        if na > 1 or ns > 0:
            return False

        if shape is None:
            return all(nonchanging(i) for i in index)
        else:
            return all(nonchanging(i, s) for i, s in zip(index, shape))
    elif isinstance(index, slice):
        return (
            index == slice(None)
            or index == slice(0, None)
            or index == slice(None, None, 1)
        )
    # slice(None) for missing dimensions in tuple index
    elif index is Ellipsis:
        return True
    elif index is np.newaxis:  # adds a dimension
        return False
    elif instance.islistgen(index):
        if shape is None:
            return False  # could be True, but we can't know when shape is not given
        else:
            if instance.isboolsequence(index):
                return all(index)
            else:
                return index == list(range(shape))
    else:
        return False


def nonchangingdims(index, ndim, axes, shape=None):
    """nonchanging for particular dimensions

    Args:
        index(index): object used in slicing (expanded)
        ndim(num): dimensions before indexings
        axes(array): dimensions for which you want to know the index
        shape(Optional(tuple)): dimension before applying index
    Returns:
        tuple
    """
    index2 = [ind for ind in index if ind is not np.newaxis]
    index2 = expand(index, ndim)
    index2 = tuple(listtools.listadvanced(index2, axes))
    if shape is not None:
        shape = tuple(listtools.listadvanced(list(shape), axes))
    b = nonchanging(index2, shape)

    axesorder, _ = axesorder_afterindexing(index, ndim)

    i = listtools.where(axesorder, lambda x: instance.islistgen(x))
    if len(i) == 1:
        i = i[0]
        if len(axesorder[i]) == 1:
            axesorder[i] = axesorder[i][0]
    try:
        b &= listtools.listadvanced(axesorder, axes) == axes
    except Exception:
        b = False
    return b


def positiveaxis(axis, ndim):
    """Positive axis

    Args:
        axis(num): dimension index
        ndim(num): number of dimensions
    Returns:
        num
    """
    if axis < 0:
        axis += ndim
    if axis < 0 or axis >= ndim:
        raise IndexError("axis out of range")
    return axis


def expand(index, ndim):
    """Expand index to ndim dimensions

    Args:
        index(index): object used in slicing
        ndim(num): number of dimensions
    Returns:
        index
    """
    if isinstance(index, tuple):
        nindex = sum([i is not np.newaxis for i in index])
        if nindex > ndim:
            raise IndexError("invalid index")

        if Ellipsis in index:
            i = index.index(Ellipsis)
            index = index[:i] + (slice(None),) * (ndim - nindex + 1) + index[i + 1 :]
        elif nindex != ndim:
            index = index + (slice(None),) * (ndim - nindex)

    elif index is Ellipsis:
        index = (slice(None),) * ndim

    elif index is np.newaxis:
        index = (np.newaxis,) + (slice(None),) * ndim

    elif ndim > 1:
        index = (index,) + (slice(None),) * (ndim - 1)

    return index


def axisindex(index, axis, ndim):
    """Find axis in expanded index

    Args:
        axis(num): dimension index
        ndim(num): number of dimensions
    Returns:
        num
    """
    axis = positiveaxis(axis, ndim)
    lst = list(np.cumsum([i is not np.newaxis for i in index]))
    return lst.index(axis + 1)


def fulldim(index, ndim):
    """Dimension after indexing (without squeezing singletons)

    Args:
        index(index):
        ndim(num): number of dimensions
    Returns:
        num
    """
    return ndim + sum([i is np.newaxis for i in index])


def replace(index, ndim, axes, rindices):
    """Replace indexing for a specified dimension

    Args:
        index(index): object used in slicing
        ndim(num): number of dimensions
        axes(list): dimension to be replaced
        rindex(list): new indexing for this dimensions
    Returns:
        index
    """
    index2 = list(expand(index, ndim))

    for axis, rindex in zip(axes, rindices):
        axis = axisindex(index2, axis, ndim)
        index2[axis] = rindex

    return tuple(index2)


def replacefull(index, ndim, axes):
    """Set the specified dimension to full range

    Args:
        index(index): object used in slicing
        ndim(num): number of dimensions
        axes(list): dimension to be set to full range
    Returns:
        index
    """
    return replace(index, ndim, axes, [slice(None)] * len(axes))


def expanddims(arr, ndim):
    """Expand array dimensions

    Args:
        arr(array):
        ndim(num): new dimensions
    Returns:
        array
    """
    n = ndim - arr.ndim
    if n > 0:
        arr = arr[(Ellipsis,) + (np.newaxis,) * n]
    return arr


class op_index(object):
    def __init__(self, index=None):
        self.index = index

    def __call__(self, data):
        if self.index is not None:
            data = data[self.index]
        return data

    def __str__(self):
        return "indexing: {}".format(self.index)


class op_transpose(object):
    def __init__(self, axes=None):
        self.axes = axes

    def __call__(self, data):
        if self.axes is not None:
            data = np.transpose(data, self.axes)

        return data

    def __str__(self):
        return "transpose: {}".format(self.axes)


class operators(object):
    def __init__(self):
        self.ops = []

    def append(self, op):
        if instance.islistgen(op):
            self.ops += op
        elif isinstance(op, operators):
            self.ops += op.ops
        else:
            if len(self.ops) > 1:
                if isinstance(self.ops[-1], op_transpose) and isinstance(
                    op, op_transpose
                ):
                    self.ops[-1].axes = listtools.listadvanced_int(
                        self.ops[-1].axes, op.axes
                    )
                else:
                    self.ops.append(op)
            else:
                self.ops.append(op)

    def __str__(self):
        return "operators:\n {}".format("\n ".join([str(o) for o in self.ops]))

    def __call__(self, data):
        for o in self.ops:
            data = o(data)
        return data


def move_axis(ndim, i1, i2):
    """Transpose axes to move axis i1 to i2

    Args:
        ndim(num):
        i1(num):
        i2(num):

    Returns:
        list
    """
    if i1 == i2:
        return None
    axes = list(range(ndim))
    axes.insert(i2, axes.pop(i1))
    return axes


def prepare_getitem(index, ndim, axis, shapefull):
    """
    Args:
        index(index): object used in slicing
        ndim(num): number of dimensions before indexing
        axis(num): generator axis
        shapefull(tuple): full stack dimension (without applying index)
    Returns:
        tuple
    """

    indexret = expand(index, ndim)

    # Remove all np.newaxis and ops2 restore operation
    index, ops2 = extract_newaxis(indexret)

    # Extract axis and ops1 restore operation
    indexrest, shaperest, indexaxis, shapeaxis, axis, ops1 = extract_axis(
        index, axis, shapefull
    )

    # Order of operations:
    #   1. np.stack([data[indexrest] ...],axis=axis)
    #   2. ops1
    #   3. ops2
    ops1.append(ops2)

    return indexret, indexrest, shaperest, indexaxis, shapeaxis, axis, ops1


def empty_afterindexing(index, ndim, shapefull=None):
    """Allocate array based on index

    Args:
        index(index): object used in slicing
        ndim(num): stack dimensions (without applying index)
        shapefull(Optional(tuple)): full stack dimension (without applying index)

    Returns:
        array
    """
    shape = shape_afterindexing(shapefull, index, ndim=ndim)
    return np.empty(shape)


def decompose_listindexing(index, nlist):
    """
    Args:
        index(index)
        nlist: resulting length of list indexing

    Returns:
        list(index)
    """

    # Replace bool with int array
    index2 = list(index)
    for i, ind in enumerate(index2):
        if instance.isboolsequence(index2):
            index2[i] = np.where(index2)

    indexlist = [None] * nlist
    nindex = len(index)
    tpl = [0] * nindex
    for j in range(nlist):
        for i, ind in enumerate(index):
            if instance.islistgen(ind):
                tpl[i] = ind[j]
            else:
                tpl[i] = ind
        indexlist[j] = tuple(tpl)

    return indexlist


def getitem(generator, arglist, index, ndim, shapefull=None, axis=0):
    """Similar to np.stack(...,axis=...)[index] except that the data
       still needs to be generate by looping over the axis dimension.
       Unnecessary data generation is avoided by inspecting index.

       >>> data = np.stack(map(generator,args),axis=axis)
       >>> shapefull = data.shape
       >>> ndim = data.ndim
       >>> data2 = getitem(generator,args,index,ndim,shapefull=shapefull,axis=axis)
       >>> np.testing.assert_array_equal(data[index],data2)

    Args:
        generator: data generator
        arglist(list): generator arguments
        index(index): object used in slicing
        ndim(num): stack dimensions (before applying index)
        shapefull(Optional(tuple)): full stack dimension (before applying index)
        axis(Optional(num)): axis corresponding to arglist

    Returns:
        array
    """

    # Better to slice before stacking, even when generator does not read from disk:
    # python -mtimeit -s'import numpy as np' 'n=50' 'np.dstack([np.zeros((n,n))[0:n//2,0:n//2]]*n)'
    # python -mtimeit -s'import numpy as np' 'n=50' 'np.dstack([np.zeros((n,n))]*n)[0:n//2,0:n//2,:]'

    # Split axis indexing from indexing of other dimensions
    index, indexrest, shaperest, indexaxis, shapeaxis, axis, ops = prepare_getitem(
        index, ndim, axis, shapefull
    )

    # Apply indexaxis on axis dimension
    bargsingleton = False
    if not nonchanging(indexaxis, shape=shapeaxis):
        if isadvanced(indexaxis):
            # No advanced indexing for list so use comprehension
            arglist = listtools.listadvanced(arglist, indexaxis)
        elif isinstance(indexaxis, numbers.Number):
            # No singleton dimension squeezing
            arglist = [arglist[indexaxis]]
            bargsingleton = True
        else:
            arglist = arglist[indexaxis]

    # Advanced indexing along axis needs to be handled separately
    if isadvanced(indexaxis):
        na = lengthadvanced(indexaxis, shapeaxis)
        nb = lengthadvanced(indexrest, shaperest)
        if nb != 0 and nb != na:
            raise IndexError(
                "shape mismatch: indexing arrays could not be broadcast together with shapes {}".format(
                    index
                )
            )
        indexlist = decompose_listindexing(indexrest, na)
        data = [generator(arg)[ind] for arg, ind in zip(arglist, indexlist)]
    else:
        data = [generator(arg)[indexrest] for arg in arglist]

    # Stack
    if len(data) == 0:
        data = empty_afterindexing(index, ndim, shapefull=shapefull)
    elif len(data) == 1 and bargsingleton:
        data = ops(data[0])
    else:
        data = np.stack(data, axis=axis)
        data = ops(data)

    return data


def setitem(
    selector, arglist, index, ndim, value, method="set", shapefull=None, axis=0
):
    """Similar to np.stack(...,axis=...)[index]=value except that the data
       still needs to be generate by looping over the axis dimension.

        Args:
        selector: data selector
        arglist(list): generator arguments
        index(index): object used in slicing
        ndim(num): stack dimensions (before applying index)
        method(Optional(str)): set, add, multiply, subtract, divide, and, or
        shapefull(Optional(tuple)): full stack dimension (without applying index)
        axis(Optional(num)): axis corresponding to arglist

    Returns:
        array
    """

    # print '-----------------------'
    # print 'index:',index
    # print 'arg axis:',axis

    # Split axis indexing from indexing of other dimensions
    index, indexrest, shaperest, indexaxis, shapeaxis, axis, ops = prepare_getitem(
        index, ndim, axis, shapefull
    )

    # Apply indexaxis on axis dimension
    if not nonchanging(indexaxis, shape=shapeaxis):
        if isadvanced(indexaxis):
            # No advanced indexing for list so use comprehension
            arglist = listtools.listadvanced(arglist, indexaxis)
        elif isinstance(indexaxis, numbers.Number):
            # No singleton dimension squeezing
            arglist = [arglist[indexaxis]]
        else:
            arglist = arglist[indexaxis]

    # Advanced indexing along axis needs to be handled separately
    shapeafter = shape_afterindexing(shapefull, index, ndim=ndim)

    if isadvanced(indexaxis):
        na = lengthadvanced(indexaxis, shapeaxis)
        nb = lengthadvanced(indexrest, shaperest)
        if nb != 0 and nb != na:
            raise IndexError(
                "shape mismatch: indexing arrays could not be broadcast together with shapes {}".format(
                    index
                )
            )
        indexlist = decompose_listindexing(indexrest, na)
    else:
        indexlist = [indexrest] * len(arglist)

    # print 'indexaxis:',indexaxis
    # print 'indexrest:',indexrest
    # print 'axis:',axis

    valueindex = not instance.isscalar(value) and axis is not None
    if valueindex:
        valueindex = [slice(None)] * len(shapeafter)
    else:
        valuesel = ops(value)

    for i, (arg, ind) in enumerate(zip(arglist, indexlist)):
        if valueindex:
            valueindex[axis] = i
            valuesel = ops(value[tuple(valueindex)])
        if method == "set":
            selector(arg)[ind] = valuesel
        elif method == "add":
            selector(arg)[ind] += valuesel
        elif method == "sub":
            selector(arg)[ind] -= valuesel
        elif method == "mul":
            selector(arg)[ind] *= valuesel
        elif method == "div":
            selector(arg)[ind] /= valuesel
        elif method == "and":
            selector(arg)[ind] &= valuesel
        elif method == "or":
            selector(arg)[ind] |= valuesel
        else:
            raise ValueError("Method {} unknown".format(method))

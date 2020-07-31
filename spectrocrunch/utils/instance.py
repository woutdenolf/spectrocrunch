# -*- coding: utf-8 -*-

import ast
import collections

try:
    import collections.abc as collections_abc
except ImportError:
    try:
        from collections import collections_abc
    except ImportError:
        collections_abc = collections
import numbers
import math
import numpy as np
import uncertainties.core
from ..patch.pint import ureg

try:
    from array import array
except ImportError:
    pass

# Builtin string types:
#                  | Python 2       |  Python 3
# -------------------------------------------------
# <byte string>     | type:str       | class:bytes
#                  | literal: ""    | literal: b""
# -------------------------------------------------
# <unicode string>  | type:unicode   | class:str
#                  | literal: u""   | literal: ""
# -------------------------------------------------
#
# <byte string>: array of bytes
# <unicode string>: array of code points
#
# code point: an integer which is given meaning by the Unicode standard
# encoding/codec (ascii, utf-8): mapping between bytes and code points
# code unit: number of bits in a code point
# grapheme: one or more code points that represent a single element
#           of the writing system
#           e.g. ä is a grapheme which, depending on the endocing,
#           is represented by one or two code points
# glyph: image to represent a grapheme or part thereof
#
# Conversion:
# encode: <unicode string> --(encoding)--> <byte string>
# decode: <byte string>    --(encoding)--> <unicode string>
#
# Implicit decoding/encoding in Python 2:
# default = sys.getdefaultencoding()
# str.encode(encoding)     -> str.decode(default).encode(encoding)
# unicode.decode(encoding) -> unicode.encode(default).decode(encoding)
#
# Builtin sequence types:
# ------------------------
#   list (mutable sequence)
#   tuple (immutable sequence)
#   bytes (immutable sequence of bytes)
#   bytearray (mutable sequence of bytes)
#   array (mutable sequence of basic C types:
#          integers of 1/2/4/8 bytes
#          floats of 4/8 bytes)
#   memoryview (byte view of sequence with buffer interface)
# Python 2:
#   str (==bytes)
#   unicode (immutable sequence code points)
#   xrange (sequence generator)
#   buffer (~=memoryview)
# Python 3:
#   str (==unicode)
#   range (==xrange)


try:
    unicode
except NameError:
    unicode = str


# Integer -> unicode
try:
    unichr
except NameError:
    unichr = chr


try:
    xrange
except NameError:
    xrange = range


listgentypes = list, xrange


def isstring(x):
    """
    Immutable sequence of code points or bytes
    """
    return isinstance(x, (bytes, unicode)) or isnpstring(x)


def asunicode(x):
    """
    Args:
        s(bytes or unicode)
    Returns:
        bytes
    Raises:
        UnicodeDecodeError: bytes with extended-ASCII characters
    """
    if isinstance(x, bytes):
        return x.decode("utf-8")
    return x


def asbytes(s):
    """
    Args:
        s(bytes or unicode)
    Returns:
        bytes
    Raises:
        UnicodeEncodeError: unicode with non-UTF8 characters
    """
    if isinstance(s, unicode):
        return s.encode("utf-8")
    return s


def issequence(x):
    """
    Sequence (mutable, immutable, generator) except for strings and bytearray's
    """
    return isinstance(x, collections_abc.Sequence) and not isinstance(
        x, (str, bytes, unicode, bytearray)
    )


def isset(x):
    return isinstance(x, collections_abc.Set)


def isiterable(x):
    """Excluding string types
    """
    if isqscalar(x):
        return False
    if isinstance(x, collections_abc.Iterable):
        return not isstring(x)
    return isqarray(x)


def isorderedmapping(x):
    return isinstance(x, collections.OrderedDict)


def ismapping(x):
    return isinstance(x, collections_abc.Mapping)


def iscallable(x):
    return isinstance(x, collections_abc.Callable)


def isarray(x):
    if issequence(x):
        return True
    types = collections_abc.Set, array, bytearray
    if isinstance(x, types):
        return True
    return isnparray(x) or isqarray(x)


def isstringarray(x):
    if not isarray(x):
        return False
    if isarray0(x):
        return isstring(x.item())
    else:
        return all(isstring(y) for y in x)


def isarray0(x):
    """numpy array with dim==0
    """
    if isnparray(x):
        return x.ndim == 0
    if isqarray(x):
        return x.magnitude.ndim == 0
    return False


def isarraynot0(x):
    """isarray, excluding numpy array with dim==0
    """
    if isarray(x):
        return not isarray0(x)
    return False


def isnparray(x):
    return isinstance(x, np.ndarray)


def isnpscalar(x):
    return np.issubclass_(x.__class__, np.generic)


def isnumpy(x):
    return isnparray(x) or isnpscalar(x)


def isnpstring(x):
    return np.issubclass_(x.__class__, np.character)


def isqarray(x):
    if isquantity(x):
        if isnparray(x.magnitude):
            return True
    return False


def isuarray(x):
    if isarray(x):
        if isarray0(x):
            return isuscalar(asscalar(x))
        else:
            return any(isuscalar(z) or isuarray(z) for z in x)
    return False


def isquantity(x):
    return isinstance(x, ureg.Quantity)


def isqscalar(x):
    if isquantity(x):
        return isscalar(x.magnitude)
    return False


def isuscalar(x):
    return isinstance(x, uncertainties.UFloat)


def isscalar(x):
    return isnumber(x) or isnpscalar(x) or isqscalar(x) or isuscalar(x)


def isnumber(x):
    return isinstance(x, numbers.Number) or isnpnumber(x)


def isnpnumber(x):
    return np.issubclass_(x.__class__, np.number)


def isnpinteger(x):
    return np.issubclass_(x.__class__, np.integer)


def isnumeric(x):
    return isnumber(x) and not isinstance(x, bool)


def isinteger(x):
    return isinstance(x, numbers.Integral) or isnpinteger(x)


def dtype_is_integer(dtype):
    return isinteger(np.array([0], dtype)[0])


def islistgen(x):
    return isinstance(x, listgentypes)


def isboolsequence(lst):
    try:
        return all(isinstance(i, bool) for i in lst) and len(lst) > 0
    except:
        return False


def asscalar(x):
    try:
        x = x.item()
    except:
        pass
    return x


class _toarray(object):
    restore = {
        "array": lambda x: x,
        "scalar": lambda x: x[0],
        "array0": lambda x: np.array(x[0]),
    }

    def __call__(self, x):
        if isquantity(x):
            # Quantity(array/scalar) -> ndarray[Quantity(scalar),...]
            u = x.units
            x, frestore = self(x.magnitude)

            def qscal(y):
                return ureg.Quantity(y, u)

            x = np.vectorize(qscal, otypes=[object])(x)

            def qrestore(y):
                return ureg.Quantity(frestore(y), u)

            return x, qrestore
        elif isarray(x):
            # Special case: numpy 0-d array
            if isnparray(x):
                if x.ndim == 0:
                    return x[np.newaxis], self.restore["array0"]
            # Create number array (possibly objects needed)
            try:
                x = np.asarray(x)
            except ValueError:
                x = np.asarray(x, dtype=object)
            return x, self.restore["array"]
        elif isscalar(x):
            return np.asarray([x]), self.restore["scalar"]
        else:
            return np.asarray([x], dtype=object), self.restore["scalar"]


asarrayf = _toarray()


def asarray(x):
    return asarrayf(x)[0]


def aslist(x):
    return asarray(x).tolist()


def asnumber(x):
    if isnumber(x):
        return x
    try:
        x = ast.literal_eval(x)
    except:
        pass
    if isnumber(x):
        return x
    else:
        try:
            return float(x)
        except:
            return np.nan


def arrayit(x):
    if isarray(x):
        return x
    else:
        return [x]

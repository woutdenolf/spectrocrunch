# -*- coding: utf-8 -*-

import random
import collections
import warnings
import numpy as np
from uncertainties import ufloat
from future.utils import with_metaclass
from abc import ABCMeta, abstractmethod, abstractproperty
from ..patch.pint import ureg

MAXDEPTH = 4
MAXITEMS = 5
NTRIES = 1000
TYPES = "numpy", "uncertainties", "pint"


try:
    unicode
except NameError:
    unicode = str


try:
    unichr
except NameError:
    unichr = chr


def try_unique_range(name, seq):
    if NTRIES:
        for i in range(NTRIES):
            yield i
        else:
            raise RuntimeError(
                "Could not find unique {} in {} tries:\n {}".format(name, NTRIES, seq)
            )
    else:
        i = 0
        while True:
            yield i
            i += 1


class DataBase(with_metaclass(ABCMeta)):
    def __init__(self, **kwargs):
        self.generate()

    @abstractmethod
    def generate(self):
        pass

    def __eq__(self, other):
        if not self.isinstance(other):
            return False
        return self.baseinstance.data == other.baseinstance.data

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return self.baseinstance.__class__.__name__

    @property
    def baseinstance(self):
        return self

    def isinstance(self, other, class_=None):
        if class_ is None:
            class_ = self.baseinstance.__class__
        return isinstance(other.baseinstance, class_)


# ====== Native misc ======


class RandomNone(DataBase):
    def generate(self):
        self.data = None


# ====== Native numbers ======


def numequal(a, b):
    # Just to make sure nan == nan
    #
    # Prevent magnitude to get converted to ndarray by np.isclose
    # Unit is always dimensionless
    try:
        a = a.magnitude
    except AttributeError:
        pass
    try:
        b = b.magnitude
    except AttributeError:
        pass
    try:
        if a == b:
            return True
    except:
        pass
    try:
        return np.isclose(a, b, rtol=0, atol=0, equal_nan=True)
    except TypeError:
        return False


class NumberBase(DataBase):
    def __eq__(self, other):
        if not self.isinstance(other, NumberBase):
            return False
        a, b = self.data, other.data
        return numequal(a, b)


class RandomBool(NumberBase):
    def generate(self):
        self.data = bool(random.getrandbits(1))


class RandomInt(NumberBase):
    def generate(self):
        self.data = random.randint(-100, 100)


def randomfloat(finite=False):
    if finite:
        return random.random()
    else:
        return random.choice([random.random()] * 5 + [float("nan"), float("inf")])


class RandomFloat(NumberBase):
    def generate(self):
        self.data = randomfloat()


class RandomComplex(NumberBase):
    def generate(self):
        self.data = complex(randomfloat(), randomfloat())


# ====== Native strings ======


class StringBase(DataBase):

    ALPHABET = None
    CLASS = None

    def __init__(self, maxitems=MAXITEMS, **kwargs):
        self.n = random.randint(0, maxitems)
        super(StringBase, self).__init__(maxitems=MAXITEMS, **kwargs)

    def generate(self):
        self.data = self._join(random.choice(self.ALPHABET) for _ in range(self.n))

    @classmethod
    @abstractmethod
    def _join(cls):
        pass

    def __eq__(self, other):
        if not self.isinstance(other, StringBase):
            return False
        a, b = self.data, other.data
        abytes = isinstance(a, bytes)
        bbytes = isinstance(b, bytes)
        if not (abytes ^ bbytes):
            if not abytes:
                a = a.encode("utf-8")
            if not bbytes:
                b = b.encode("utf-8")
        return a == b

    def __repr__(self):
        s = super(StringBase, self).__repr__()
        return "{}_{}".format(s, self.n)


def alphabet_bytes(encoding="ascii"):
    include_ranges = [(0x0021, 0x0021), (0x0023, 0x0026), (0x0028, 0x007E)]
    if encoding != "ascii":
        include_ranges += [(0x00A1, 0x00AC), (0x00AE, 0x00FF)]
    alphabet = [
        unichr(code_point).encode(encoding)
        for current_range in include_ranges
        for code_point in range(current_range[0], current_range[1] + 1)
    ]
    return alphabet


class RandomBytes(StringBase):

    CLASS = bytes

    @classmethod
    def _join(cls, iterable):
        return cls.CLASS(b"".join(iterable))


class RandomBytesAscii(RandomBytes):

    ALPHABET = alphabet_bytes(encoding="ascii")


class RandomBytesExt(RandomBytes):

    ALPHABET = alphabet_bytes(encoding="latin1")


def alphabet_unicode():
    include_ranges = [
        (0x0021, 0x0021),
        (0x0023, 0x0026),
        (0x0028, 0x007E),
        (0x00A1, 0x00AC),
        (0x00AE, 0x00FF),
        (0x0100, 0x017F),
        (0x0180, 0x024F),
        (0x2C60, 0x2C7F),
        (0x16A0, 0x16F0),
        (0x0370, 0x0377),
        (0x037A, 0x037E),
        (0x0384, 0x038A),
        (0x038C, 0x038C),
    ]
    alphabet = [
        unichr(code_point)
        for current_range in include_ranges
        for code_point in range(current_range[0], current_range[1] + 1)
    ]
    return alphabet


class RandomUnicode(StringBase):

    ALPHABET = alphabet_unicode()
    CLASS = unicode

    @classmethod
    def _join(cls, iterable):
        return cls.CLASS(u"".join(iterable))


# ====== Native sequences ======


def init_sequence(
    seq_types,
    unique=False,
    _depth=0,
    nitems=None,
    maxdepth=MAXDEPTH,
    maxitems=MAXITEMS,
    **kwargs
):
    """
    Generate a sequence with random length, random types
    Args:
        seq_types: possible item types
        unique: sequence should have unique values
        maxdepth: maximal sequence depth
        maxitems: maximal number of items
    """
    if nitems is None:
        if _depth == maxdepth:
            n = 0
        else:
            n = random.randint(0, maxitems)
    else:
        n = nitems
    kwargs = kwargs.copy()
    kwargs["_depth"] = _depth + 1
    kwargs["maxdepth"] = maxdepth
    kwargs["maxitems"] = maxitems
    if unique:
        # Make sure the sequence could (but not always) have unique values
        seq = []
        reprs = collections.Counter()
        reprmax = {"bool": 2, "emptystring": 1}
        reprmap = {
            "RandomBool": "bool",
            "RandomNumpyBool": "bool",
            "RandomBytesAscii_0": "emptystring",
            "RandomBytesExt_0": "emptystring",
            "RandomUnicode_0": "emptystring",
        }
        for _ in try_unique_range("sequence", seq):
            if len(seq) == n:
                break
            item = random.choice(seq_types)(**kwargs)
            repritem = repr(item)
            repritem = reprmap.get(repritem, repritem)
            nrepr = reprmax.get(repritem, None)
            if item not in seq and not reprs[repritem] == nrepr:
                seq.append(item)
                reprs[repritem] += 1
    else:
        seq = [random.choice(seq_types)(**kwargs) for _ in range(n)]
    return seq


def generate_sequence(seq, unique=False):
    if unique:
        for _ in try_unique_range("generate_sequence", seq):
            for item in seq:
                item.generate()
            useq = []
            for item in seq:
                if item in useq:
                    break
                useq.append(item)
            else:
                break
    else:
        for item in seq:
            item.generate()


class SequenceBase(DataBase):

    CLASS = None

    def __init__(self, **kwargs):
        if not hasattr(self, "_values"):
            self._values = init_sequence(self.seq_types(), unique=False, **kwargs)
        super(SequenceBase, self).__init__(**kwargs)

    @classmethod
    def seq_types(cls):
        return (RandomData,)

    def __eq__(self, other):
        if not self.isinstance(other):
            return False
        return self._values == other._values

    @property
    def data(self):
        return self.CLASS(v.data for v in self._values)

    def generate(self):
        generate_sequence(self._values, unique=False)


class HashableSequenceBase(SequenceBase):
    @classmethod
    def seq_types(cls):
        return (RandomHashable,)


class RandomTuple(SequenceBase):

    CLASS = tuple


class RandomHashableTuple(HashableSequenceBase):

    CLASS = tuple


class RandomList(SequenceBase):

    CLASS = list


# ====== Native sets ======


def sort_equal(a, b):
    if len(a) != len(b):
        return False
    for v in a:
        if v not in b:
            return False
    return True


class SetBase(HashableSequenceBase):
    def __init__(self, **kwargs):
        self._values = init_sequence(self.seq_types(), unique=True, **kwargs)
        super(SetBase, self).__init__(**kwargs)

    def __eq__(self, other):
        if not self.isinstance(other):
            return False
        return sort_equal(self._values, other._values)

    @property
    def data(self):
        values = [v.data for v in self._values]
        if values:
            random.shuffle(values)
        return self.CLASS(values)

    def generate(self):
        generate_sequence(self._values, unique=True)


class RandomSet(SetBase):

    CLASS = set


class RandomFrozenSet(SetBase):

    CLASS = frozenset


# ====== Native mappings ======


class OrderedMappingBase(DataBase):

    CLASS = None

    def __init__(self, **kwargs):
        keys = init_sequence((RandomHashable,), unique=True, **kwargs)
        self._values = init_sequence(
            (RandomData,), unique=False, nitems=len(keys), **kwargs
        )
        self._keys = keys

    def __eq__(self, other):
        if not self.isinstance(other):
            return False
        return self._keys == other._keys and self._values == other._values

    def generate(self):
        generate_sequence(self._keys, unique=True)
        generate_sequence(self._values, unique=False)

    def data_items(self):
        keys = [k.data for k in self._keys]
        values = [k.data for k in self._values]
        return zip(keys, values)

    @property
    def data(self):
        return self.CLASS(self.data_items())


class UnorderedMappingBase(OrderedMappingBase):
    def __eq__(self, other):
        if not self.isinstance(other):
            return False
        if not sort_equal(self._keys, other._keys):
            return False
        for k, v in zip(self._keys, self._values):
            i = other._keys.index(k)
            if v != other._values[i]:
                return False
        return True

    def data_items(self):
        keys = [k.data for k in self._keys]
        values = [k.data for k in self._values]
        items = list(zip(keys, values))
        random.shuffle(items)
        return items


class RandomDict(UnorderedMappingBase):

    CLASS = dict


class RandomOrderedDict(OrderedMappingBase):

    CLASS = collections.OrderedDict


# ====== Numpy ======


class NumpyScalarBase(NumberBase):
    def __init__(self, **kwargs):
        self._class = random.choice(self.NPTYPES)
        super(NumpyScalarBase, self).__init__(**kwargs)

    @classmethod
    @abstractmethod
    def datagen(cls):
        pass

    def generate(self):
        self.data = self._class(*self.datagen())


class RandomNumpyBool(NumpyScalarBase):

    NPTYPES = (np.bool,)

    @classmethod
    def datagen(cls):
        return (random.getrandbits(1),)


class RandomNumpyInt(NumpyScalarBase):

    NPTYPES = (
        np.byte,
        np.short,
        np.int,
        np.longlong,
        np.int8,
        np.int16,
        np.int32,
        np.int64,
    )

    @classmethod
    def datagen(cls):
        return (random.randint(-100, 100),)


class RandomNumpyUInt(NumpyScalarBase):

    NPTYPES = (
        np.ubyte,
        np.ushort,
        np.uint,
        np.ulonglong,
        np.uint8,
        np.uint16,
        np.uint32,
        np.uint64,
    )

    @classmethod
    def datagen(cls):
        return (random.randint(0, 100),)


class RandomNumpyFloat(NumpyScalarBase):

    NPTYPES = (np.single, np.double, np.float, np.float16, np.float32, np.float64)

    @classmethod
    def datagen(cls):
        return (randomfloat(),)


class RandomNumpyComplex(NumpyScalarBase):

    NPTYPES = (np.complex,)

    @classmethod
    def datagen(cls):
        return randomfloat(), randomfloat()


class RandomNumpyArray(SequenceBase):

    CLASS = np.array

    @classmethod
    def seq_types(cls):
        return (RandomNumber,)

    @property
    def data(self):
        return self.CLASS([v.data for v in self._values])


class RandomNumpyArray0(RandomNumpyArray):
    def __init__(self, **kwargs):
        super(RandomNumpyArray0, self).__init__(nitems=1, **kwargs)

    @property
    def data(self):
        return self.CLASS(self._values[0].data)


# ====== Uncertainties ======


class RandomUNumber(NumpyScalarBase):

    NPTYPES = (ufloat,)

    @classmethod
    def datagen(cls):
        return randomfloat(finite=True), randomfloat(finite=True)


# ====== Pint ======


class RandomPintNumber(NumberBase):

    UNIT = ureg.dimensionless

    def __init__(self, **kwargs):
        self._magnitude = RandomPintMagnitude(**kwargs)
        super(RandomPintNumber, self).__init__(**kwargs)

    def generate(self):
        self._magnitude.generate()

    @property
    def data(self):
        return ureg.Quantity(self._magnitude.data, self.UNIT)


class RandomPintArray(SequenceBase):

    CLASS = ureg.Quantity
    UNIT = ureg.dimensionless

    @classmethod
    def seq_types(cls):
        return (RandomPintMagnitude,)

    @property
    def data(self):
        magnitudes = [v.data for v in self._values]
        return self.CLASS(magnitudes, self.UNIT)


class RandomPintArray0(RandomPintArray):
    def __init__(self, **kwargs):
        super(RandomPintArray0, self).__init__(nitems=1, **kwargs)

    @property
    def data(self):
        magnitude = self._values[0].data
        return self.CLASS(np.array(magnitude), self.UNIT)


# ====== Choices ======


class ChoiceBase(DataBase):
    def __init__(self, **kwargs):
        self.data = random.choice(self.choices(**kwargs))(**kwargs)
        super(ChoiceBase, self).__init__(**kwargs)

    @classmethod
    @abstractmethod
    def choices(cls, **kwargs):
        pass

    def __eq__(self, other):
        return self.baseinstance == other.baseinstance

    @property
    def data(self):
        return self._choice.data

    @data.setter
    def data(self, value):
        self._choice = value

    def generate(self):
        self._choice.generate()

    def __repr__(self):
        return repr(self._choice)

    @property
    def baseinstance(self):
        return self._choice.baseinstance


class RandomNativeNumber(ChoiceBase):
    @classmethod
    def choices(cls, **kwargs):
        return (RandomBool, RandomInt, RandomFloat, RandomComplex)


class RandomNumpyNumber(ChoiceBase):
    @classmethod
    def choices(cls, **kwargs):
        return (RandomNumpyBool, RandomNumpyInt, RandomNumpyFloat, RandomNumpyComplex)


class RandomHashableNumber(ChoiceBase):
    @classmethod
    def choices(cls, types=TYPES, **kwargs):
        ret = [RandomNativeNumber]
        if "numpy" in types:
            ret.append(RandomNumpyNumber)
        return tuple(ret)


class RandomNumber(ChoiceBase):
    @classmethod
    def choices(cls, types=TYPES, **kwargs):
        ret = [RandomNativeNumber]
        if "numpy" in types:
            ret.append(RandomNumpyNumber)
        if "uncertainties" in types:
            ret.append(RandomUNumber)
        if "pint" in types:
            ret.append(RandomPintNumber)
        return tuple(ret)


class RandomPintMagnitude(ChoiceBase):
    @classmethod
    def choices(cls, types=TYPES, **kwargs):
        ret = RandomInt, RandomFloat
        if "numpy" in types:
            ret += RandomNumpyInt, RandomNumpyFloat
        if "uncertainties" in types:
            ret += (RandomUNumber,)
        return ret


class RandomString(ChoiceBase):
    @classmethod
    def choices(cls, **kwargs):
        return RandomBytesAscii, RandomBytesExt, RandomUnicode


class RandomAtom(ChoiceBase):
    @classmethod
    def choices(cls, **kwargs):
        return (
            (RandomNone,)
            + RandomNumber.choices(**kwargs)
            + RandomString.choices(**kwargs)
        )


class RandomHashableAtom(ChoiceBase):
    @classmethod
    def choices(cls, **kwargs):
        return (
            (RandomNone,)
            + RandomHashableNumber.choices(**kwargs)
            + RandomString.choices(**kwargs)
        )


class RandomMutableSequence(ChoiceBase):
    @classmethod
    def choices(cls, **kwargs):
        return RandomList, RandomSet


class RandomImmutableSequence(ChoiceBase):
    @classmethod
    def choices(cls, **kwargs):
        return RandomTuple, RandomFrozenSet


class RandomNumpySequence(ChoiceBase):
    @classmethod
    def choices(cls, **kwargs):
        return RandomNumpyArray, RandomNumpyArray0


class RandomPintSequence(ChoiceBase):
    @classmethod
    def choices(cls, **kwargs):
        return RandomPintArray, RandomPintArray0


class RandomSequence(ChoiceBase):
    @classmethod
    def choices(cls, **kwargs):
        return (
            RandomMutableSequence,
            RandomImmutableSequence,
            RandomNumpySequence,
            RandomPintSequence,
        )


class RandomHashableSequence(ChoiceBase):
    @classmethod
    def choices(cls, **kwargs):
        return RandomHashableTuple, RandomFrozenSet


class RandomHashable(ChoiceBase):
    @classmethod
    def choices(cls, **kwargs):
        return RandomHashableAtom.choices(**kwargs) + RandomHashableSequence.choices(
            **kwargs
        )


class RandomMapping(ChoiceBase):
    @classmethod
    def choices(cls, **kwargs):
        return RandomDict, RandomOrderedDict


class RandomData(ChoiceBase):
    @classmethod
    def choices(cls, **kwargs):
        return RandomAtom, RandomSequence, RandomMapping


def factory(**kwargs):
    return RandomData(**kwargs)

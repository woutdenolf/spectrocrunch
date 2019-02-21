import random
import string
import collections
import numpy as np
from future.utils import with_metaclass
from abc import ABCMeta, abstractmethod, abstractproperty




MAXDEPTH = 4
MAXITEMS = 5


class DataBase(with_metaclass(ABCMeta)):

    def __init__(self, **kwargs):
        self.generate(**kwargs)

    @abstractmethod
    def generate(self, **kwargs):
        pass

    def __eq__(self, other):
        if isinstance(other, DataBase):
            other = other.data
        return self.data == other

    def __ne__(self, other):
        return not self.__eq__(other)

    def __str__(self):
        return types.unicode(self.data)


# ====== Not iterables ======


class RandomNone(DataBase):

    def generate(self, **kwargs):
        self.data = None


class RandomBool(DataBase):

    def generate(self, **kwargs):
        self.data = bool(random.getrandbits(1))


class RandomInt(DataBase):

    def generate(self, **kwargs):
        self.data = random.randint(-100, 100)


class RandomFloat(DataBase):

    def generate(self, **kwargs):
        self.data = random.random()


class StringBase(DataBase):

    ALPHABET = None
    CLASS = None

    def generate(self, maxitems=MAXITEMS, **kwargs):
        n = random.randint(0, maxitems)
        self.data = self._join(random.choice(self.ALPHABET) for _ in range(n))

    @classmethod
    @abstractmethod
    def _join(cls):
        pass


def alphabet_bytes(encoding='ascii'):
    try:
        get_char = unichr
    except NameError:
        get_char = chr
    include_ranges = [
        (0x0021, 0x0021),
        (0x0023, 0x0026),
        (0x0028, 0x007E),
    ]
    if encoding != 'ascii':
        include_ranges += [
            (0x00A1, 0x00AC),
            (0x00AE, 0x00FF)
        ]
    alphabet = [get_char(code_point).encode(encoding)
        for current_range in include_ranges
        for code_point in range(current_range[0], current_range[1] + 1)
    ]
    return alphabet


class RandomBytes(StringBase):

    CLASS = bytes

    @classmethod
    def _join(cls, iterable):
        return cls.CLASS(b''.join(iterable))


class RandomAscii(RandomBytes):

    ALPHABET = alphabet_bytes(encoding='ascii')


class RandomLatin1(RandomBytes):

    ALPHABET = alphabet_bytes(encoding='latin1')


def alphabet_unicode():
    try:
        get_char = unichr
    except NameError:
        get_char = chr
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
    alphabet = [get_char(code_point)
        for current_range in include_ranges
        for code_point in range(current_range[0], current_range[1] + 1)
    ]
    return alphabet


class RandomUnicode(StringBase):

    ALPHABET = alphabet_unicode()
    try:
        CLASS = unicode
    except NameError:
        CLASS = str

    @classmethod
    def _join(cls, iterable):
        return cls.CLASS(u''.join(iterable))


# ====== Sequences and sets ======


def init_sequence(seq_types, unique=False, _depth=0, nitems=None,
                  maxdepth=MAXDEPTH, maxitems=MAXITEMS):
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
    kwargs = {'_depth': _depth+1,
              'maxdepth': maxdepth,
              'maxitems': maxitems}
    if unique:
        # Make sure the sequence could (but not always) have unique values
        _seq = []
        seq = []
        while len(seq) != n:
            value = random.choice(seq_types)(**kwargs)
            _value = value.data
            if _value not in _seq:
                _seq.append(_value)
                seq.append(value)
    else:
        seq = [random.choice(seq_types)(**kwargs) for _ in range(n)]
    return seq


class SequenceBase(DataBase):

    CLASS = None

    def __init__(self, **kwargs):
        self._values = init_sequence(self.seq_types(),
                                     unique=False,
                                     **kwargs)
        super(SequenceBase, self).__init__()

    @classmethod
    def seq_types(cls):
        return RandomData,

    @property
    def data(self):
        return self.CLASS(v.data for v in self._values)
        
    def generate(self, **kwargs):
        for item in self._values:
            item.generate(**kwargs)


class HashableSequenceBase(SequenceBase):

    @classmethod
    def seq_types(cls):
        return RandomHashable,


class SetBase(DataBase):

    def __init__(self, **kwargs):
        self._values = init_sequence(self.seq_types(),
                                     unique=True,
                                     **kwargs)
        super(SetBase, self).__init__()

    @classmethod
    def seq_types(cls):
        return RandomHashable,

    @property
    def data(self):
        values = [v.data for v in self._values]
        if values:
            random.shuffle(values)
        return self.CLASS(values)
        
    def generate(self, **kwargs):
        # Values need to be unique
        while True:
            for item in self._values:
                item.generate(**kwargs)
            if len(self.data) == len(self._values):
                break


class RandomTuple(SequenceBase):

    CLASS = tuple


class RandomHashableTuple(HashableSequenceBase):

    CLASS = tuple


class RandomList(SequenceBase):

    CLASS = list


class RandomSet(SetBase):

    CLASS = set


class RandomFrozenSet(SetBase):

    CLASS = frozenset


class RandomNumpy(SequenceBase):

    CLASS = np.array

    @classmethod
    def seq_types(cls):
        return RandomNumber,


# ====== Mappings ======


class OrderedMappingBase(DataBase):

    CLASS = None

    def __init__(self, **kwargs):
        keys = init_sequence((RandomHashable,),
                             unique=True,
                             **kwargs)
        self._values = init_sequence((RandomData,),
                                     unique=False,
                                     nitems=len(keys),
                                     **kwargs)
        self._keys = keys

    def generate(self, **kwargs):
        # Keys need to be unique
        while True:
            n = len(self._keys)
            for item in self._keys:
                item.generate(**kwargs)
            keys = [k.data for k in self._keys]
            if len(set(keys)) == n:
                break
        for item in self._values:
            item.generate(**kwargs)

    def data_items(self):
        keys = [k.data for k in self._keys]
        values = [k.data for k in self._values]
        return zip(keys, values)

    @property
    def data(self):
        return self.CLASS(self.data_items())


class UnorderedMappingBase(OrderedMappingBase):

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


# ====== Choices======


class ChoiceBase(DataBase):

    def __init__(self, **kwargs):
        self.data = random.choice(self.choices())(**kwargs)
        super(ChoiceBase, self).__init__()

    @classmethod
    @abstractmethod
    def choices(cls):
        pass

    @property
    def data(self):
        return self._choice.data

    @data.setter
    def data(self, value):
        self._choice = value

    def generate(self, **kwargs):
        self._choice.generate(**kwargs)

    def __repr__(self):
        return repr(self._choice)


class RandomNumber(ChoiceBase):

    @classmethod
    def choices(cls):
        return RandomBool, RandomInt, RandomFloat


class RandomString(ChoiceBase):

    @classmethod
    def choices(cls):
        return RandomAscii, RandomLatin1, RandomUnicode


class RandomAtom(ChoiceBase):

    @classmethod
    def choices(cls):
        return (RandomNone,) + RandomNumber.choices() +\
                RandomString.choices()


class RandomSequence(ChoiceBase):

    @classmethod
    def choices(cls):
        return RandomTuple, RandomList, RandomSet, RandomFrozenSet


class RandomHashableSequence(ChoiceBase):

    @classmethod
    def choices(cls):
        return RandomHashableTuple, RandomFrozenSet


class RandomHashable(ChoiceBase):

    @classmethod
    def choices(cls):
        return RandomAtom.choices() + RandomHashableSequence.choices()


class RandomMapping(ChoiceBase):

    @classmethod
    def choices(cls):
        return RandomDict, RandomOrderedDict


class RandomData(ChoiceBase):

    @classmethod
    def choices(cls):
        return RandomAtom, RandomSequence, RandomMapping


def factory(maxdepth=MAXDEPTH, maxitems=MAXITEMS):
    return RandomData(maxdepth=maxdepth, maxitems=maxitems)

# -*- coding: utf-8 -*-

import re
import logging
from ..utils import hashing


logger = logging.getLogger(__name__)


def locate_number(name):
    start, end, num = 0, 0, '0'
    for m in re.finditer('[0-9]+', name):
        start, end, num = m.start(), m.end(), m.group()
    return start, end, int(num), start == end


def naming_format(name):
    start, end, num, nonumber = locate_number(name)
    if nonumber:
        fmt = name+'.{:d}'
    else:
        fmt = name[:start]+'{{:0{}d}}'.format(end-start)+name[end:]
    return fmt, num, nonumber


def regex(name):
    name = re.escape(name)
    start, end, num, nonumber = locate_number(name)
    if nonumber:
        pattern = '^' + name + r'(\.[0-9]+)?$'
    else:
        pattern = '^' + name[:start]+'[0-9]+' + name[end:] + '$'
    return pattern


def match(name):
    return re.compile(regex(name)).match


def increment(name, add=1):
    """
    name -> name.1
    name.1 -> name.2
    name0001 -> name0002
    name0001_0001 -> name0001_0002
    """
    fmt, num, _ = naming_format(name)
    return fmt.format(num+add)


def prepare(name):
    """
    name -> name.1
    name0001 -> name0001
    """
    fmt, num, nonumber = naming_format(name)
    if nonumber:
        num = 1
    return fmt.format(num)


class Name(object):

    def __init__(self, name):
        fmt, num, nonumber = naming_format(name)
        if nonumber:
            num = 1
        self.fmt = fmt
        self.num = num
        self._orgnum = num
        self.locked = False

    def reset(self):
        self.num = self._orgnum
    
    def lock(self):
        self.locked = True

    def unlock(self):
        self.locked = True

    def __repr__(self):
        return repr(str(self))

    def __str__(self):
        return self.fmt.format(self.num)

    def __iadd__(self, x):
        if self.locked:
            logger.warning('{} is locked and therefore not incremented'
                .format(repr(self)))
        else:
            self.num += x
        return self

    def __add__(self, x):
        return self.__class__(self.fmt.format(self.num+x))

    @property
    def matchfunc(self):
        return match(str(str))


def calc_checksum(dependencies, parameters, alreadyhash=False):
    if dependencies:
        hashes = [prev.checksum for prev in dependencies]
    else:
        hashes = []
    if alreadyhash:
        h = parameters
    else:
        h = hashing.calchash(parameters)
    hashes.append(h)
    return hashing.mergehash(*hashes)

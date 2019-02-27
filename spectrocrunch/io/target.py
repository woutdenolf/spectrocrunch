# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
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

import re

from ..utils import hashing


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
        pattern = '^'+name+'(\.[0-9]+)?$'
    else:
        pattern = '^'+name[:start]+'[0-9]+'+name[end:]+'$'
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

    def reset(self):
        self.num = self._orgnum

    def __repr__(self):
        return self.fmt.format(self.num)

    def __str__(self):
        return repr(self)

    def __iadd__(self, x):
        self.num += x
        return self

    def __add__(self, x):
        return self.__class__(self.fmt.format(self.num+1))

    @property
    def matchfunc(self):
        return match(str(str))


def calc_checksum(dependencies, confighash):
    if dependencies:
        hashes = [prev.checksum for prev in dependencies]
    else:
        hashes = []
    hashes.append(confighash)
    return hashing.mergehash(*hashes)

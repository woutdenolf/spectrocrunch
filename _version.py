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
"""
Module for version handling (adopted from silx)
* version = "1.2.3" or "1.2.3-beta4"
* version_info = named tuple (1,2,3,"beta",4)
* hexversion: 0x010203B4
* strictversion = "1.2.3b4
"""

MAJOR = 0
MINOR = 0
MICRO = 1
RELEV = "dev"
SERIAL = 1

RELEASE_LEVEL_VALUE = {"dev": 0,
                       "alpha": 10,
                       "beta": 11,
                       "rc": 12,
                       "final":15}

from collections import namedtuple
_version_info = namedtuple("version_info", ["major", "minor", "micro", "releaselevel", "serial"])

version_info = _version_info(MAJOR, MINOR, MICRO, RELEV, SERIAL)

strictversion = version = "%d.%d.%d" % version_info[:3]

if version_info.releaselevel != "final":
    version += "-%s%s" % version_info[-2:]
    prerel = "a" if RELEASE_LEVEL_VALUE.get(version_info[3], 0) < 10 else "b"
    if prerel not in "ab":
        prerel = "a"
    strictversion += prerel + str(version_info[-1])

hexversion = version_info[4]
hexversion |= RELEASE_LEVEL_VALUE.get(version_info[3], 0) * 1 << 4
hexversion |= version_info[2] * 1 << 8
hexversion |= version_info[1] * 1 << 16
hexversion |= version_info[0] * 1 << 24
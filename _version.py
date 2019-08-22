# -*- coding: utf-8 -*-

from collections import namedtuple
"""
Module for version handling (adopted from silx)
* version = "1.2.3" or "1.2.3-beta4"
* version_info = named tuple (1,2,3,"beta",4)
* hexversion: 0x010203B4
* strictversion = "1.2.3b4"
"""

MAJOR = 0
MINOR = 0
MICRO = 3  # <=15
RELEV = "dev"
SERIAL = 1  # <=15

RELEASE_LEVEL_VALUE = {"dev": 0,
                       "alpha": 10,
                       "beta": 11,
                       "rc": 12,
                       "final": 15}

_version_info = namedtuple(
    "version_info", ["major", "minor", "micro", "releaselevel", "serial"])

version_info = _version_info(MAJOR, MINOR, MICRO, RELEV, SERIAL)

strictversion = version = "%d.%d.%d" % version_info[:3]

if version_info.releaselevel != "final":
    version += "-%s%s" % version_info[-2:]
    prerel = "a" if RELEASE_LEVEL_VALUE.get(version_info[3], 0) <= 10 else "b"
    strictversion += prerel + str(version_info[-1])

hexversion = version_info[4]
hexversion |= RELEASE_LEVEL_VALUE.get(version_info[3], 0) * 1 << 4
hexversion |= version_info[2] * 1 << 8
hexversion |= version_info[1] * 1 << 16
hexversion |= version_info[0] * 1 << 24

if __name__ == "__main__":
    print(version)

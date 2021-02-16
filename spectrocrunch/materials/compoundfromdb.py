# -*- coding: utf-8 -*-

from . import compoundfromname
from . import compoundfromnist
from . import compoundfromformula
from . import compoundsearch

import collections
import itertools

db = collections.OrderedDict(
    [
        ("name", compoundfromname),
        ("nist", compoundfromnist),
        ("formula", compoundfromformula),
    ]
)


def getnames():
    return list(
        itertools.chain(compoundfromname.getnames(), compoundfromnist.getnames())
    )


def factory(name, **kwargs):
    for dbname, mod in db.items():
        match = mod.search(name)
        if match:
            return mod.factory(match[0], **kwargs)
    return None


def search(name):
    return compoundsearch.search(getnames(), name)


def compoundfromdb(name):
    return factory(name)

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

from collections import OrderedDict

from . import instance

import future.utils

# https://blog.ionelmc.ro/2015/02/09/understanding-python-metaclasses/
#
# Instance creation:
#   metaclass.__call__()
#       instance = class.__new__(self,)
#       class.__init__(instance,...)
#       return instance
#
# Class creation:
#   metametaclass.__call__()
#       class = metaclass.__new__(self,)
#       metaclass.__init__(class,...)
#       return class


def clsfactory(cls, name):
    """
    Get class from name/alias

    Args:
        name(str): name of class that needs to be created

    Returns:
        class
    """
    if name in cls.clsregistry:
        return cls.clsregistry[name]
    elif name in cls.aliasregistry:
        return cls.aliasregistry[name]
    else:
        raise RuntimeError("Class {} is not known:\n registered classes: {}\n aliases: {}".format(
            name, cls.clsnames(), cls.clsaliases()))


def factory(cls, name, *args, **kwargs):
    """
    Get class instance from class name/alias

    Args:
        name(str): name of class that needs to be created

    Returns:
        class instance
    """
    return cls.clsfactory(name)(*args, **kwargs)


def register(cls, regcls, name):
    """
    Register class name and aliases

    Args:
        regcls(class): class to be registered
        name(str): name of the class to be registered
    """
    lname = name.lower()
    cls.clsregistry[name] = regcls
    if lname != name:
        cls.aliasregistry[lname] = regcls
    if hasattr(regcls, "aliases"):
        for alias in regcls.aliases:
            lalias = alias.lower()
            cls.aliasregistry[alias] = regcls
            if lalias != alias:
                cls.aliasregistry[lalias] = regcls


def clsnames(cls):
    """Registered class names
    """
    return list(cls.clsregistry.keys())


def clsaliases(cls):
    """Registered class aliases
    """
    return list(cls.aliasregistry.keys())


def clsallnames(cls):
    """Registered class names+aliases
    """
    return cls.clsnames()+cls.clsaliases()


class FactoryMeta(type):
    """
    Metaclass used to register all classes inheriting from FactoryMeta 
    """

    def __new__(self, name, bases, attr):
        cls = super(FactoryMeta, self).__new__(self, name, bases, attr)

        if not hasattr(cls, "register"):
            cls.clsregistry = OrderedDict()
            cls.aliasregistry = OrderedDict()
            cls.clsnames = classmethod(clsnames)
            cls.clsallnames = classmethod(clsallnames)
            cls.clsaliases = classmethod(clsaliases)
            cls.clsfactory = classmethod(clsfactory)
            cls.factory = classmethod(factory)
            cls.register = classmethod(register)

        return cls

    def __init__(cls, name, bases, attr):
        cls.register(cls, name)
        for b in bases:
            if hasattr(b, "register"):
                b.register(cls, name)
        super(FactoryMeta, cls).__init__(name, bases, attr)


def with_metaclass(bases=None):
    if bases is None:
        return future.utils.with_metaclass(FactoryMeta)
    else:
        if not instance.isarray(bases):
            bases = (bases,)
        return future.utils.with_metaclass(FactoryMeta, *bases)

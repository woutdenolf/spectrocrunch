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

from six import with_metaclass

class FactoryMeta(type):
    """
    Metaclass used to register all lens classes inheriting from FactoryBase 
    """
    def __init__(cls, name, bases, dct):
        if name!="FactoryBase":
            lname = name.lower()
            cls.registry[name] = cls
            if lname!=name:
                cls.registry2[lname] = cls
            if hasattr(cls, "aliases"):
                for alias in cls.aliases:
                    lalias = alias.lower()
                    cls.registry2[alias] = cls
                    if lalias!=alias:
                        cls.registry2[lalias] = cls
        super(FactoryMeta, cls).__init__(name, bases, dct)

class FactoryBase(with_metaclass(FactoryMeta, object)):
    """
    Class representing a factory class
    """

    @classmethod
    def factory(cls, name, **kwargs):
        """
        Args:
            name(str): class name

        Returns:
            object: class instance
        """
        if name in cls.registry:
            return cls.registry[name](**kwargs)
        elif name in cls.registry2:
            return cls.registry2[name](**kwargs)
        else:
            raise RuntimeError("Class {} is not one of the registered classes: {}".format(name, cls.registry.keys()+cls.registry2.keys()))

def ClassFactory(name):
    """
    Example of a factory function (not useful by itself)
    """
    def __init__(self, **kwargs):
        FactoryBase.__init__(self, **kwargs)
    newclass = type(name, (FactoryBase,),{"__init__": __init__})
    return newclass


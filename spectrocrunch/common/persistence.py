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

import json
import importlib
from . import instance

class SerializedGenerator(dict):

    def __init__(self,module=None,generator=None,args=tuple(),kwargs=dict()):
        if not instance.isstring(module):
            raise RuntimeError("Serialized generator must specify a module")
        if not instance.isstring(generator):
            raise RuntimeError("Serialized generator must specify a generator")
        self["name"] = "spectrocrunch.serialized"
        self["module"] = module
        self["generator"] = generator
        self["args"] = args
        self["kwargs"] = kwargs


def deserialize(data,top=True):
    if instance.isstring(data) and top:
        with open(data,'r') as f:
            data = json.load(f)

    if isinstance(data,dict):
        if data.get("name")=="spectrocrunch.serialized":
            m = importlib.import_module(data["module"])
            
            g = data["generator"]
            if not hasattr(m,g):
                raise RuntimeError("The module {} does not have the attribute {} anymore".format(m,g))
            g = getattr(m,g)
            
            args = tuple(deserialize(arg,top=False) for arg in data["args"])
            kwargs = {k:deserialize(v,top=False) for k,v in data["kwargs"].items()}
            return g(*args,**kwargs)

        for k in data:
            data[k] = deserialize(data[k],top=False)
            
    return data


def serialize(data):

    if hasattr(data,"serialize"):
        data = data.serialize()
    elif hasattr(data,"units"):
        from .units import serialize as func
        data = func(data)

    return data 
    

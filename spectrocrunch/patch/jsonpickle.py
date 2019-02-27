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

from __future__ import absolute_import
import jsonpickle
from jsonpickle import Pickler, Unpickler
import jsonpickle.ext.numpy as jsonpickle_numpy
from jsonpickle.handlers import BaseHandler, register
from .pint import ureg


jsonpickle.set_preferred_backend('json')
jsonpickle.set_encoder_options('json', sort_keys=True)
# jsonpickle.set_decoder_options('json', ...)


try:
    unicode
except NameError:
    unicode = str


def dumps(data):
    return jsonpickle.encode(data, keys=True)


def loads(data):
    return jsonpickle.decode(data, keys=True)


def flatten(data, **kwargs):
    """
    Args:
        data(Any):
        \**kwargs: see Pickler
    Returns:
        dict: state of data
    """
    return Pickler(**kwargs).flatten(data, reset=True)


def restore(state, **kwargs):
    """
    Args:
        state(dict):
        \**kwargs: see Pickler
    Returns:
        Any: data restored from state
    """
    return Unpickler.restore(state, reset=True, **kwargs)


class QuantityHandler(BaseHandler):

    def flatten(self, obj, state):
        tpl = obj.to_tuple()
        tpl = self.context.flatten(tpl, reset=False)
        state['tpl'] = tpl
        return state

    def restore(self, state):
        tpl = self.context.restore(state['tpl'], reset=False)
        return ureg.Quantity.from_tuple(tpl)


class UnitHandler(BaseHandler):

    def flatten(self, obj, state):
        state['units'] = str(obj)
        return state

    def restore(self, state):
        return ureg.Unit(state['units'])


jsonpickle_numpy.register_handlers()
register(ureg.Quantity, QuantityHandler, base=True)
register(ureg.Unit, UnitHandler, base=True)

# -*- coding: utf-8 -*-

from __future__ import absolute_import
import jsonpickle
import jsonpickle.ext.numpy as jsonpickle_numpy
from jsonpickle.handlers import BaseHandler, register
from .pint import ureg


jsonpickle.set_preferred_backend('json')
jsonpickle.set_encoder_options('json', sort_keys=True)


try:
    unicode
except NameError:
    unicode = str


def dumps(data, **kwargs):
    jsonpickle.set_encoder_options('json', **kwargs)
    return jsonpickle.encode(data, keys=False)


def loads(data, **kwargs):
    jsonpickle.set_decoder_options('json', **kwargs)
    return jsonpickle.decode(data, keys=False)


def flatten(data, **kwargs):
    """
    Args:
        data(Any):
        \**kwargs: see Pickler
    Returns:
        dict: state of data
    """
    return jsonpickle.Pickler(**kwargs).flatten(data, reset=True)


def restore(state, **kwargs):
    """
    Args:
        state(dict):
        \**kwargs: see Pickler
    Returns:
        Any: data restored from state
    """
    return jsonpickle.Unpickler.restore(state, reset=True, **kwargs)


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

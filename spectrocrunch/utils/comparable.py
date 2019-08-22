# -*- coding: utf-8 -*-

import operator


class Comparable(object):

    @property
    def _repr(self):
        """Unique representation of an instance
        """
        return "{}{}".format(type(self).__name__, id(self))

    def _cmpkey(self, other):
        return self._repr

    def _sortkey(self, other):
        return self._cmpkey

    def __repr__(self):
        return self._repr

    def __str__(self):
        return repr(self)

    def encode(self, *args, **kwargs):
        return str(self).encode(*args, **kwargs)

    def _compareop(self, other, op, key):
        a = getattr(self, key)(other)
        try:
            b = getattr(other, key)(self)
        except:
            return op(a, other)
        else:
            return op(a, b)

    def _sort(self, other, op):
        return self._compareop(other, op, '_sortkey')

    def _compare(self, other, op):
        return self._compareop(other, op, '_cmpkey')

    def __lt__(self, other):
        return self._sort(other, operator.lt)

    def __le__(self, other):
        return self._sort(other, operator.le)

    def __ge__(self, other):
        return self._sort(other, operator.ge)

    def __gt__(self, other):
        return self._sort(other, operator.gt)

    def __eq__(self, other):
        return self._compare(other, operator.eq)

    def __ne__(self, other):
        return self._compare(other, operator.ne)

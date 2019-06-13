# -*- coding: utf-8 -*-
#
#   Copyright (C) 2016 European Synchrotron Radiation Facility, Grenoble, France
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

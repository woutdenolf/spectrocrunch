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

import operator


class Function(object):
    """Lazy function with standard logical and arithmetic operators.
    All lazy functions will be evaluated with the same arguments.
    """

    def __init__(self, func, name=None):
        self._func = func
        self._name = name

    def __str__(self):
        name = self._name
        if not name:
            name = self._func.__name__
        return name

    def __call__(self, *args, **kwargs):
        return self._func(*args, **kwargs)

    @classmethod
    def _binary_op(cls, a, b, op, symb1, symb2=None):
        lazya = isinstance(a, cls)
        lazyb = isinstance(b, cls)
        if lazya and lazyb:
            f = lambda *args, **kwargs: op(a(*
                                             args, **kwargs), b(*args, **kwargs))
        elif lazya:
            f = lambda *args, **kwargs: op(a(*args, **kwargs), b)
        elif lazyb:
            f = lambda *args, **kwargs: op(a, b(*args, **kwargs))
        else:
            f = lambda *args, **kwargs: op(a, b)

        if symb2:
            return cls(f, name="{}{},{}{}".format(symb1, a, b, symb2))
        else:
            return cls(f, name="({}{}{})".format(a, symb1, b))

    @classmethod
    def _unitary_op(cls, a, op, symb1, symb2=None):
        if isinstance(a, cls):
            f = lambda *args, **kwargs: op(a(*args, **kwargs))
        else:
            f = lambda *args, **kwargs: op(a)

        if symb2:
            return cls(f, name="{}{}{}".format(symb1, a, symb2))
        else:
            return cls(f, name="({}{})".format(symb1, a))

    def __mul__(self, rother):
        return self._binary_op(self, rother, operator.mul, '*')

    def __div__(self, rother):
        return self._binary_op(self, rother, operator.truediv, '/')

    def __truediv__(self, rother):
        return self._binary_op(self, rother, operator.truediv, '/')

    def __floordiv__(self, rother):
        return self._binary_op(self, rother, operator.floordiv, '//')

    def __mod__(self, rother):
        return self._binary_op(self, rother, operator.mod, '%')

    def __divmod__(self, rother):
        return self._binary_op(self, rother, operator.divmod, 'divmod(', ')')

    def __add__(self, rother):
        return self._binary_op(self, rother, operator.add, '+')

    def __sub__(self, rother):
        return self._binary_op(self, rother, operator.sub, '-')

    def __pow__(self, rother):
        return self._binary_op(self, rother, operator.pow, '^')

    def __rmul__(self, lother):
        return self._binary_op(lother, self, operator.mul, '*')

    def __rdiv__(self, lother):
        return self._binary_op(lother, self, operator.truediv, '/')

    def __rtruediv__(self, lother):
        return self._binary_op(lother, self, operator.truediv, '/')

    def __rfloordiv__(self, lother):
        return self._binary_op(lother, self, operator.floordiv, '//')

    def __rmod__(self, lother):
        return self._binary_op(lother, self, operator.mod, '%')

    def __rdivmod__(self, lother):
        return self._binary_op(lother, self, operator.mod, 'divmod(', ')')

    def __radd__(self, lother):
        return self._binary_op(lother, self, operator.add, '+')

    def __rsub__(self, lother):
        return self._binary_op(lother, self, operator.sub, '-')

    def __rpow__(self, lother):
        return self._binary_op(lother, self, operator.pow, '^')

    def __lt__(self, rother):
        return self._binary_op(self, rother, operator.lt, '<')

    def __le__(self, rother):
        return self._binary_op(self, rother, operator.le, '<=')

    def __gt__(self, rother):
        return self._binary_op(self, rother, operator.gt, '>')

    def __ge__(self, rother):
        return self._binary_op(self, rother, operator.ge, '>=')

    def __eq__(self, rother):
        return self._binary_op(self, rother, operator.eq, '==')

    def __ne__(self, rother):
        return self._binary_op(self, rother, operator.ne, '!=')

    def __and__(self, rother):
        return self._binary_op(self, rother, operator.and_, ' and ')

    def __or__(self, rother):
        return self._binary_op(self, rother, operator.or_, ' or ')

    def __not__(self):
        return self._unitary_op(self, operator.not_, 'not ')

    def __neg__(self):
        return self._unitary_op(self, operator.neg, '-')

    def __pos__(self):
        return self._unitary_op(self, operator.pos, '+')

    def __abs__(self):
        return self._unitary_op(self, operator.abs, 'abs(', ')')

    # def __float__(self):
    #    return self._unitary_op(self,float,'float(',')')

    # def __long__(self):
    #    return self._unitary_op(self,long,'long(',')')

    # def __int__(self):
    #    return self._unitary_op(self,int,'int(',')')

    # def __bool__(self):
    #    return self._unitary_op(self,operator.truth,'==True')

    # def __nonzero__(self):
    #    return self._binary_op(self,0,operator.ne,"!=")

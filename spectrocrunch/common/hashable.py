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

from .comparable import Comparable

class Hashable(Comparable):
    """A derived class with methods _cmpkey() and _stringrepr() can be used
       as dictionary keys. Additionally, the string representations can still
       be used as keys.
    """

    def _compare(self, other, method):
        if isinstance(other,str):
            return method(str(self), other)
        else:
            return method(self._cmpkey(), other._cmpkey())

    def _sort(self, other, method):
        if isinstance(other,str):
            return method(str(self), other)
        else:
            return method(self._sortkey(), other._sortkey())
            
    def _cmpkey(self):
        return self.__str__()

    def _sortkey(self):
        return self.__str__()
        
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        try:
            return self._stringrepr()
        except (AttributeError):
            return str(id(self))
            
    def __hash__(self):
        return hash(self.__str__())

    




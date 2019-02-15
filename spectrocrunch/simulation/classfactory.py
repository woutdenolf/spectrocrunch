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

from ..utils.classfactory import FactoryMeta
from ..utils import instance

import future.utils
import numpy as np


class SimulClass(object):

    @staticmethod
    def propagate_broadcast(N, *args):
        """
        Args:
            N(num|array): incomming number of photons
            args(tuple(num|array)): energy related variables

        Returns:
            unumpy.uarray: len(energy) x len(N)
        """
        if instance.isarray(N) or instance.isarray(args[0]):
            nN = np.asarray(N).shape
            if len(nN) == 2:
                nN = nN[1]
            else:
                nN = int(np.product(nN))
            nenergy = np.asarray(args[0]).size
            N = np.broadcast_to(N, [nenergy, nN])
            args = tuple(np.broadcast_to(arg, [nN, nenergy]).T for arg in args)
        return (N,)+args


def with_metaclass(bases=None):
    if bases is None:
        return future.utils.with_metaclass(FactoryMeta, SimulClass)
    else:
        if not instance.isarray(bases):
            bases = (bases,)
        return future.utils.with_metaclass(FactoryMeta, SimulClass, *bases)

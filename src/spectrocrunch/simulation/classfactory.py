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
                nN = int(np.prod(nN))
            nenergy = np.asarray(args[0]).size
            N = np.broadcast_to(N, [nenergy, nN])
            args = tuple(np.broadcast_to(arg, [nN, nenergy]).T for arg in args)
        return (N,) + args


def with_metaclass(bases=None):
    if bases is None:
        return future.utils.with_metaclass(FactoryMeta, SimulClass)
    else:
        if not instance.isarray(bases):
            bases = (bases,)
        return future.utils.with_metaclass(FactoryMeta, SimulClass, *bases)

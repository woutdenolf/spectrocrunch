# -*- coding: utf-8 -*-

from copy import copy, deepcopy


class Copyable(object):

    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result

    def __deepcopy__(self, memo):
        selfid = id(self)
        if selfid in memo:
            return memo[selfid]
        else:
            cls = self.__class__
            result = cls.__new__(cls)
            memo[id(self)] = result
            for k, v in self.__dict__.items():
                setattr(result, k, deepcopy(v, memo))
            return result

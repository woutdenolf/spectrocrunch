# -*- coding: utf-8 -*-

class Enum(set):

    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError

    def __call__(self, name):
        return self.__getattr__(name)

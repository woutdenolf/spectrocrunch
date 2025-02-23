from .comparable import Comparable


class Hashable(Comparable):
    @property
    def _repr(self):
        return "{}{}".format(type(self).__name__, id(self))

    def __hash__(self):
        return hash(self._repr)


class CompHashable(Hashable, Comparable):
    @property
    def _repr(self):
        return "{}{}".format(type(self).__name__, id(self))

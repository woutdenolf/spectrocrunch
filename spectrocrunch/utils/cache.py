# -*- coding: utf-8 -*-

from collections import OrderedDict
import contextlib
import functools


class LimitedSizeDict(OrderedDict):
    
    def __init__(self, *args, **kwds):
        self.size_limit = kwds.pop("size_limit", None)
        OrderedDict.__init__(self, *args, **kwds)
        self._check_size_limit()

    def __setitem__(self, key, value):
        OrderedDict.__setitem__(self, key, value)
        self._check_size_limit()

    def _check_size_limit(self):
        if self.size_limit is not None:
            while len(self) > self.size_limit:
                self.popitem(last=False)


class Cache(object):

    def __init__(self, force=False):
        self._store = {}
        self.force = force

    def cachegenerator(self, subject):
        method = "_cache_{}".format(subject)
        try:
            generator = getattr(self, method)
        except AttributeError:
            generator = None
        return generator

    def _generate_cache(self, subject, *args, **kwargs):
        generator = self.cachegenerator(subject)
        if generator is None:
            raise RuntimeError(
                "Cache subject \"{}\" does not exist".format(subject))
        return generator(*args, **kwargs)

    @contextlib.contextmanager
    def cachectx(self, subject, *args, **kwargs):
        """Context manager to cache a subject with a particular function
        """

        # Already cashed?
        cache = subject not in self._store

        # Cache
        if cache:
            self._store[subject] = self._generate_cache(
                subject, *args, **kwargs)

        # Use cache
        yield

        # Keep cache?
        if cache:
            self._store.pop(subject, None)

    def getcache(self, subject):
        """Get cache belonging to a subject
        """
        if subject in self._store:
            return self._store[subject]
        else:
            if self.force:
                if self.cachegenerator(subject) is not None:
                    raise RuntimeError(
                        "Cache \"{}\" can only be used in context \" self.cachectx(\"{}\") \"".format(subject, subject))
                else:
                    raise AttributeError("'{}' object has no attribute '{}'".format(
                        self.__class__.__name__, subject))
            else:
                # generate the context
                return self._generate_cache(subject)


def withcache(subject):
    """Decorator to cache a subject
    """
    def cache(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            try:
                ctx = getattr(self, "cachectx")
            except AttributeError:
                raise RuntimeError(
                    "Class \"{}\" must be derived from \"Cache\"".format(self.__class__.__name__))
            with ctx(subject):
                return func(self, *args, **kwargs)
        return wrapper
    return cache

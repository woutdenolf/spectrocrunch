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

import os
import contextlib
from abc import ABCMeta, abstractproperty, abstractmethod
from future.utils import with_metaclass
import operator
import logging
import functools
import re

from . import utils
from ..utils import instance

logger = logging.getLogger(__name__)


class FileSystemException(Exception):
    """
    Base class for generic file system exceptions.
    """
    pass


class Missing(FileSystemException):
    """
    Raised when a file system operation can't be performed because
    a file or directory does not exist when it is required to.
    """
    pass


class AlreadyExists(FileSystemException):
    """
    Raised when a file system operation can't be performed because
    a file or directory exists but is required to not exist.
    """
    pass


class MissingParentDirectory(FileSystemException):
    """
    Raised when a parent directory doesn't exist.
    (Imagine mkdir without -p)
    """
    pass


class NotADirectory(FileSystemException):
    """
    Raised when a file system operation can't be performed because
    an expected directory is actually a file.
    """
    pass


class NotAFile(FileSystemException):
    """
    Raised when a file system operation can't be performed because
    an expected file is actually a directory.
    """
    pass


class DirectoryIsNotEmpty(FileSystemException):
    """
    Raised when a file system operation can't be performed because
    a directory is not empty.
    """
    pass


class File(with_metaclass(ABCMeta, object)):
    """Proxy to file
    """

    def __init__(self, **kwargs):
        self._handle = None
        self._onclose_callbacks = []
        super(File, self).__init__()

    @contextlib.contextmanager
    def open(self, **openparams):
        if self._handle is None:
            self._openparams_defaults(openparams)
            with self._closectx():
                with self._fopen(**openparams) as self._handle:
                    try:
                        yield self._handle
                    finally:
                        self._handle = None
        else:
            yield self._handle

    def _openparams_defaults(self, openparams):
        for k, v in self.openparams.items():
            if k not in openparams:
                openparams[k] = v

    @contextlib.contextmanager
    @abstractmethod
    def _fopen(self, **openparams):
        pass

    @abstractproperty
    def openparams(self):
        pass

    @contextlib.contextmanager
    def _closectx(self):
        try:
            yield
        finally:
            if self._onclose_callbacks:
                funcs, self._onclose_callbacks = self._onclose_callbacks, []
                for func, args, kwargs in funcs:
                    func(*args, **kwargs)

    def isopen(self):
        # only checks local handle (could be opened by others)
        return bool(self._handle)

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __repr__(self):
        pass


def onclose(func):
    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        if self._handle:
            self._onclose_callbacks.append((func, (self,)+args, kwargs))
        else:
            return func(self, *args, **kwargs)
    return wrapper


class Path(File):
    """Proxy to file system path
    """

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, value):
        self._path = self.norm(value)

    @property
    def location(self):
        if self.device:
            return '{}{}{}'.format(self.device, self.devsep, self.path)
        else:
            return self.path

    @property
    def uri(self):
        return self.location

    @property
    def device(self):
        return ''

    @property
    def devsep(self):
        return ''

    def devsplit(self, path):
        return path.split(self.devsep)

    def _split_path(self, path, device=None):
        """Split device and path.

        .. code-block:: 

            path = '/tmp/test.h5:/entry/subentry:name'
            device = None
            _split_path(path,device)==('/tmp/test.h5','/entry/subentry:name')

            path = '/tmp/test.h5:/entry/subentry:name'
            device = '/tmp/test.h5'
            _split_path(path,device)==('/tmp/test.h5','/entry/subentry:name')

            path = '/tmp/test.h5:/entry/subentry:name'
            device = '/tmp/other.h5'
            _split_path(path,device)==('/tmp/other.h5','/tmp/test.h5:/entry/subentry:name')

        Args:
            path(str):
            device(Optional(str)):

        Returns:
            device(str):
            path(str):
        """
        path = str(path)
        if self.devsep:
            # Tricky because a path could also have self.devsep characters
            tmp = self.devsplit(path)
            n = len(tmp)
            if n == 1:
                # /tmp/test.h5
                devicepath, localpath = tmp[0], ''
            elif n == 2:
                # /tmp/test.h5:/entry/subentry
                devicepath, localpath = tmp
            else:
                # /tmp/test.h5:/entry/subentry:name
                devicepath, localpath = tmp[0], self.devsep.join(tmp[1:])

            if device is not None:
                strdevice = str(device)
                if devicepath == strdevice:
                    devicepath = device
                else:
                    if path.startswith(strdevice):
                        devicepath = device
                        localpath = path[len(strdevice):]
                        if localpath.startswith(self.devsep):
                            localpath = localpath[len(self.devsep):]
                    elif not localpath:
                        devicepath = device
                        localpath = path
        else:
            devicepath, localpath = None, path

        return devicepath, localpath

    @property
    def factory_kwargs(self):
        return {}

    @property
    def factorycls(self):
        return self.__class__

    def factory(self, path):
        cls = self.factorycls
        if isinstance(path, cls):
            return path
        else:
            return cls(path, **self.factory_kwargs)

    def detach(self, **kwargs):
        return self.factorycls(self.location, **kwargs)

    def readonly(self, **kwargs):
        kwargs['mode'] = 'r'
        self.detach(**kwargs)

    def sibling(self, path):
        if isinstance(path, Path):
            return path
        _device, _path = self._split_path(path, device=self.device)
        if _device is None:
            _device = ''
        else:
            _device = str(_device)
        if str(self.device) == _device and not self.isabs(_path):
            return self.parent[_path]
        else:
            return self.factory(path)

    @staticmethod
    def _getpath(path):
        try:
            return path.path
        except AttributeError:
            return path

    def __repr__(self):
        return repr(str(self))

    def __str__(self):
        return self.location

    def __hash__(self):
        return hash(str(self))

    def __bool__(self):
        return bool(str(self))

    def __nonzero__(self):
        return self.__bool__()

    def __contains__(self, other):
        return self[other].exists

    def _binary_op(self, other, op):
        return op(str(self), str(other))

    def __eq__(self, other):
        return self._binary_op(other, operator.eq)

    def __ne__(self, other):
        return self._binary_op(other, operator.ne)

    def __le__(self, other):
        return self._binary_op(other, operator.le)

    def __ge__(self, other):
        return self._binary_op(other, operator.ge)

    def __lt__(self, other):
        return self._binary_op(other, operator.lt)

    def __gt__(self, other):
        return self._binary_op(other, operator.gt)

    def __len__(self):
        return len(str(self))

    def __add__(self, other):
        return self[other]

    def __iadd__(self, other):
        self.path += other

    def __getitem__(self, value):
        # Do not use self.location in self.join
        path = self.join(self.path, value)
        if self.device:
            path = "{}{}{}".format(self.device, self.devsep, path)
        return self.factory(path)

    #def __call__(self, *value):
    #    return self.__getitem__(value)

    def __iter__(self):
        for path in self.listdir():
            yield path

    def append(self, value):
        self += value

    @property
    def isfile(self):
        return not self.isdir

    @property
    def parent(self):
        return self['..']

    def common(self, path):
        path = self.factory(path)
        if path.device != self.device:
            return None
        lst1 = [p for p in self.iterup]
        lst2 = [p for p in path.iterup]
        ret = None
        for p1, p2 in zip(reversed(lst1), reversed(lst2)):
            if p1 == p2:
                ret = p1
            else:
                return ret
        return ret

    def _copy_move_prepare(self, dest, force=False, rename=False):
        """Called by `move` and `copy`
        """
        dest = self.sibling(dest)
        if dest.isdir and not rename:
            dest = dest[self.name]
        if dest.exists:
            if force:
                dest.remove(recursive=True)
            else:
                raise AlreadyExists(dest.location)
        d = dest.parent
        if d and not d.exists:
            d.mkdir()
        return dest

    def _move_copyrenamedel(self, dest, rename=False):
        """Instead of move, do copy(temporary name) + rename + delete
        """
        tmp = dest.parent.randomnode(prefix=dest.name)
        self.copy(tmp, force=True, follow=True)
        ret = tmp.rename(dest.path)
        self.remove(recursive=True)
        return ret

    def _mkdir_prepare(self, recursive=True, force=True):
        if self.exists:
            if not force:
                raise AlreadyExists(self.location)
            elif self.isfile:
                raise NotADirectory(self.location)
            else:
                return True
        return False

    @staticmethod
    def _link_follow(lnk):
        while lnk:
            dest = lnk.linkdest(follow=True)
            if dest:
                lnk = dest
            else:
                break
        return lnk

    @abstractproperty
    def exists(self):
        # For links: destination exists
        pass

    @abstractproperty
    def lexists(self):
        # For links: link exists
        pass

    @abstractproperty
    def isdir(self):
        pass

    @abstractproperty
    def islink(self):
        pass

    @abstractmethod
    def listdir(self):
        pass

    @abstractproperty
    def sep(self):
        pass

    def _sep_in(self, path):
        if os.sep == self.sep:
            return path
        else:
            return path.replace(self.sep, os.sep)

    def _sep_out(self, path):
        if os.sep == self.sep:
            return path
        else:
            osroot = os.path.abspath(os.sep)
            if path.startswith(osroot):
                path = os.sep + path[len(osroot):]
            return path.replace(os.sep, self.sep)

    def join(self, *args):
        args = list(map(self._sep_in, args))
        return self._sep_out(os.path.join(*args))

    def norm(self, path):
        """remove redundant separators
           remove references '.' and '..'
           make path absolute
        """
        path = self._sep_in(path)
        path = os.path.abspath(os.path.normpath(path))
        return self._sep_out(path)

    def abspath(self, path):
        return self.norm(self.join(self.parent.path, path))

    def relpath(self, path):
        """get path relative to self
        """
        path = self._sep_in(path)
        if os.path.isabs(path):
            path = os.path.relpath(path, self.path)
        return self._sep_out(path)

    def isabs(self, path):
        """absolute or relative path
        """
        path = self._sep_in(path)
        return os.path.isabs(path)

    @property
    def name(self):
        return os.path.basename(self._sep_in(self.path))

    @property
    def root(self):
        for root in self.iterup:
            pass
        return root

    @property
    def iterup(self):
        # yields self as well
        path = self
        while path != path.parent and path.parent:
            yield path
            path = path.parent
        yield path

    def tree(self, recursive=True, depth=0, files=True, _level=0):
        """
        Returns:
            2-tuple: (Path|None, list(2-tuple)|None)
        """
        if self.isdir:
            lst = []
            out = (self, lst)

            for path in sorted(self.listdir()):
                if path.isdir:
                    if not recursive and _level == depth:
                        lst.append((path, None))
                    else:
                        lst.append(path.tree(recursive=recursive, depth=depth,
                                             files=files, _level=_level+1))
                elif files:
                    lst.append((path, None))
        else:
            if files:
                out = (self, None)
            else:
                out = (None, None)
        return out

    def strtree(self, recursive=True, depth=0, files=True, stats=False):
        tree = self.tree(recursive=recursive, depth=depth, files=files)
        lines = self._str_tree(tree, stats=stats)
        return '\n'.join(lines)

    def ls(self, recursive=False, depth=0, files=True, stats=False):
        tree = self.strtree(recursive=recursive, depth=depth,
                            files=files, stats=stats)
        print(tree)

    @classmethod
    def _str_tree(cls, node, stats=False, _padding=' '):
        path, lst = node
        out = []
        if path:
            # Node name
            indent = _padding[:-1]
            if indent:
                node = path.name
            else:
                node = path.location
            nodename = '{}+-{}'.format(indent, node)

            # Add link info
            lnkdest = path.linkdest(follow=False)
            if lnkdest:
                if lnkdest.lexists:
                    sep = '-'
                else:
                    sep = '/'
                nodename = '{} -{}-> '.format(nodename, sep)

                if path.device == lnkdest.device:
                    lnkdeststr = path.parent.relpath(lnkdest.path)
                    lst = None
                else:
                    lnkdeststr = lnkdest.location
                    _padding += ' '*(len(nodename)-4)

                nodename += lnkdeststr
                if stats:
                    nodename += lnkdest._contentinfo()
            else:
                if path.exists:
                    if stats:
                        nodename += path._contentinfo()
                else:
                    nodename = nodename + ' (does not exist)'
            out.append(nodename)

            # Add stats
            if stats:
                try:
                    sts = path.stats()
                except Missing:
                    pass  # broken link
                else:
                    if lst:
                        sep = '|'
                    else:
                        sep = ' '
                    for k, v in sts.items():
                        out.append('{} {} @{} = {}'.format(
                            _padding, sep, k, v))

            # Add subnodes
            if lst:
                last = len(lst)
                _padding = _padding + ' '
                for i, node in enumerate(lst, 1):
                    out.append(_padding+'|')
                    if i == last:
                        out += cls._str_tree(node, stats=stats,
                                             _padding=_padding+' ')
                    else:
                        out += cls._str_tree(node, stats=stats,
                                             _padding=_padding+'|')
        return out

    def _contentinfo(self):
        return ''

    def node(self, follow=True):
        if follow:
            lnk = self.linkdest(follow=True)
            if lnk:
                node = lnk
            else:
                node = self
        else:
            node = self
        return node

    @abstractmethod
    def mkdir(self, recursive=True, force=True):
        pass

    def mkfile(self, **params):
        with self.open(**params):
            pass
        return self

    @abstractmethod
    def read(self):
        pass

    @abstractmethod
    def write(self, *args, **kwargs):
        pass

    @abstractmethod
    def move(self, dest, force=False, rename=False):
        """
        Atomic if the underlying filesystem has an atomic `rename`.
        Filesystem boundaries are handled with copytemp/delete/rename.
        """
        pass

    def mv(self, *args, **kwargs):
        return self.move(*args, **kwargs)

    def rename(self, dest, force=False):
        return self.move(dest, force=force, rename=True)

    @abstractmethod
    def copy(self, dest, force=False, follow=False, dereference=False):
        pass

    def cp(self, *args, **kwargs):
        return self.copy(*args, **kwargs)

    @abstractmethod
    def stats(self, follow=True):
        pass

    def get_stat(self, key, default=None, follow=True):
        try:
            return self.stats(follow=follow).get(key, default=default)
        except fs.Missing:
            return default

    @abstractmethod
    def update_stats(self, follow=True, **stats):
        pass

    @contextlib.contextmanager
    def openstats(self, follow=True):
        stats = self.stats(follow=follow)
        yield stats
        self.update_stats(**stats)

    @abstractmethod
    def remove(self, recursive=False):
        pass

    def rm(self, *args, **kwargs):
        return self.remove(*args, **kwargs)

    def rmdir(self):
        if self.isdir:
            self.remove(recursive=True)

    @abstractmethod
    def link(self, dest, soft=True):
        pass

    @abstractmethod
    def linkdest(self, follow=False):
        pass

    @property
    def linkdestname(self):
        lnkdest = self.linkdest(follow=False)
        if lnkdest:
            if self.device == lnkdest.device:
                lnkdest = self.relpath(lnkdest.path)
            else:
                lnkdest = lnkdest.location
        return lnkdest

    @contextlib.contextmanager
    def temp(self, name=None, force=False, **kwargs):
        """Context manager which creates a non-existing temporary path
        that will be removed or renamed on exit.
        """
        path = self.randomnode(**kwargs)
        try:
            yield path
        except Exception as e:
            path.remove(recursive=True)
            raise e
        finally:
            if name:
                if path.exists:
                    path.renameremove(name, force=force)
            else:
                path.remove(recursive=True)

    def renameremove(self, dest, force=False):
        """Rename or remove if the destination already exists (unless force=True)
        """
        try:
            return self.rename(dest, force=force)
        except AlreadyExists:
            self.remove(recursive=True)
            return self

    def randomnode(self, prefix='', suffix='', **kwargs):
        """Non-existing node with random name
        """
        path = self[prefix+utils.randomstring(**kwargs)+suffix]
        while path.exists:
            path = self[prefix+utils.randomstring(**kwargs)+suffix]
        return path

    def find(self, match, recursive=False, files=True, depth=0, _level=0):
        """Find files/directories

        Args:
            match(str|callable): regex when str
        """
        if instance.isstring(match):
            match = re.compile(match).match
        for path in self:
            if path.isfile and not files:
                continue
            if match(path):
                yield path
            if recursive and path.isdir and _level > depth:
                for sub in path.find(match,
                                     recursive=True,
                                     depth=depth,
                                     _level=_level+1):
                    yield sub

    def nonexisting_name(self, name):
        if name not in self:
            return name
        m = re.match(r'^(.+)\((\d+)\)$', name)
        if m:
            name = m.groups()[0]
            i = int(m.groups()[1])+1
        else:
            i = 1
        fmt = name + '({})'
        while fmt.format(i) in self:
            i += 1
        return fmt.format(i)

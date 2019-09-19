# -*- coding: utf-8 -*-

import os
import errno
import shutil
from contextlib import contextmanager
import logging
import tempfile
from datetime import datetime

from . import fs
from . import utils

logger = logging.getLogger(__name__)


class FileSystemException(fs.FileSystemException):
    """
    Base class for generic local file system exceptions.
    """
    pass


class Path(fs.Path):
    """Proxy to local path
    """

    def __init__(self, path, mode='a+', **kwargs):
        """
        Args:
            path(str):
            mode(Optional(str)):
                r : open (error when not existing), read from beginning
                r+: open (error when not existing), read/write from beginning
                w : create (truncate existing), write
                w+: create (truncate existing), read/write
                x : create (error when existing), write
                x+: create (error when existing), read/write
                a : create (append when existing), write from end
                a+: create (append when existing), read/write from end
        """
        if mode not in ('r', 'r+', 'w', 'w+', 'x', 'x+', 'a', 'a+'):
            raise ValueError('Invalid mode {}'.format(repr(mode)))
        self.openparams = kwargs
        self.openparams['mode'] = mode

        self.path = path
        super(Path, self).__init__(**kwargs)

    @property
    def mode(self):
        return self.openparams['mode']

    @property
    def factory_kwargs(self):
        return self.openparams

    @property
    def openparams(self):
        return self._openparams

    @openparams.setter
    def openparams(self, value):
        self._openparams = value

    def _openparams_defaults(self, openparams):
        super(Path, self)._openparams_defaults(openparams)
        defaultmode = self.openparams['mode']
        mode = openparams['mode']
        if defaultmode == 'r':
            # read-only
            mode = defaultmode
        elif defaultmode == 'r+' and mode not in ['r', 'r+']:
            # deny new files
            mode = defaultmode
        elif defaultmode in ['x', 'x+'] and mode not in ['w', 'w+']:
            # allow new files (error when existing)
            mode = defaultmode
        elif defaultmode in ['w', 'w+'] and mode not in ['r', 'w', 'w+']:
            # allow new files (overwrite when existing)
            mode = defaultmode
        elif defaultmode in ['a', 'a+'] and mode not in ['r', 'r+', 'a', 'a+']:
            # allow new files (append when existing)
            mode = defaultmode
        openparams['mode'] = mode

    @contextmanager
    def _fopen(self, **openparams):
        #msg = '{} ({})'.format(self.path,openparams)
        try:
            #logger.debug('Open '+msg)
            with open(self.path, **openparams) as f:
                yield f
        except IOError as err:
            if err.errno == errno.ENOENT:
                raise fs.Missing(self.location)
            elif err.errno == errno.EISDIR:
                raise fs.NotAFile(self.location)
            elif not self.isfile:
                # On windows EACCES instead of EISDIR
                raise fs.NotAFile(self.location)
            else:
                raise
        # finally:
        #    logger.debug('Close '+msg)

    @property
    def exists(self):
        return os.path.exists(self.path)

    @property
    def lexists(self):
        return os.path.lexists(self.path)

    @property
    def isdir(self):
        return os.path.isdir(self.path)

    @property
    def isfile(self):
        return os.path.isfile(self.path)

    @property
    def islink(self):
        return os.path.islink(self.path)

    @property
    def directory(self):
        return self.factory(os.path.dirname(self.path))

    def listdir(self):
        if self.isdir:
            for k in os.listdir(self.path):
                yield self.factory(self.join(self.path, k))

    @property
    def sep(self):
        return os.sep

    def mkdir(self, recursive=True, force=True):
        if self._mkdir_prepare(recursive=recursive, force=force):
            return self
        if recursive:
            try:
                os.makedirs(self.path)
            except OSError as err:
                if err.errno != errno.EEXIST:
                    raise
        else:
            parent = self.parent
            if not parent.exists:
                raise fs.MissingParentDirectory(self.location)
            os.mkdir(self.path)
        return self

    def mkfile(self, data=None, **params):
        with self.open(**params) as f:
            if data:
                f.write(data)
        return self

    def move(self, dest, force=False, rename=False):
        dest = self._copy_move_prepare(dest, force=force, rename=rename)
        try:
            os.rename(self.path, dest.path)
        except OSError as err:
            if err.errno == errno.EXDEV:
                dest = self._move_copyrenamedel(dest)
            else:
                raise
        return dest

    def copy(self, dest, force=False, follow=False, dereference=False):
        dest = self._copy_move_prepare(dest, force=force)
        if self.islink and not follow:
            # just copy the link
            dest.link(self.linkdest)
        else:
            # TODO: missing option: expand softlinks when pointing outside
            #       the tree being copied.
            if self.isdir:
                shutil.copytree(self.path, dest.path, symlinks=not dereference)
            else:
                shutil.copy(self.path, dest.path)
        return dest

    def remove(self, recursive=False):
        if self.islink:
            if recursive:
                node = self.linkdest()
                if node:
                    node.remove(recursive=True)
            os.remove(self.path)
        elif self.isdir:
            if recursive:
                shutil.rmtree(self.path)
            else:
                try:
                    os.rmdir(self.path)
                except OSError as err:
                    if err.errno == errno.ENOTEMPTY:
                        raise fs.DirectoryIsNotEmpty(self.location)
                    else:
                        raise
        elif self.exists:
            os.remove(self.path)

    def stats(self, follow=True):
        try:
            if follow:
                ret = os.stat(self.path)
            else:
                ret = os.lstat(self.path)
        except OSError as err:
            if err.errno == errno.ENOENT:
                raise fs.Missing(self.location)

        ret = {k[3:]: getattr(ret, k) for k in dir(ret) if k.startswith('st_')}

        for k, v in ret.items():
            if 'time' in k:
                ret[k] = datetime.fromtimestamp(v)

        #ret['permissions'] = oct(ret.pop('mode') & 0o777)
        return ret

    def update_stats(self, follow=True, mode=None, uid=-1, gid=-1, **ignore):
        if mode is not None:
            #mode = int(mode, 8)
            if follow:
                os.chmod(self.path, mode)
            else:
                os.lchmod(self.path, mode)
        if uid >= 0 or gid >= 0:
            if follow:
                os.chown(self.path, uid, gid)
            else:
                os.lchown(self.path, uid, gid)

    def link(self, dest, soft=True):
        lnkname = self.path
        if soft:
            dest = self.parent.relpath(self._getpath(dest))
            os.symlink(dest, lnkname)
        else:
            dest = self.sibling(dest).path
            os.link(dest, lnkname)
        return self

    def linkdest(self, follow=False):
        try:
            lnk = self.factory(self.abspath(os.readlink(self.path)))
            if follow:
                lnk = self._link_follow(lnk)
            return lnk
        except OSError:
            return None

    def read(self):
        with self.open(mode='r') as f:
            return f.read()

    def write(self, content, **openparams):
        with self.open(**openparams) as f:
            f.write(content)
        return self


@contextmanager
def temp(**kwargs):
    """
    Context manager which creates a non-existing temporary path that
    will be removed or renamed on exit.

    Args:
        **kwargs: see localfs.Path.temp
    """
    with Path(tempfile.gettempdir()).temp(**kwargs) as path:
        yield path

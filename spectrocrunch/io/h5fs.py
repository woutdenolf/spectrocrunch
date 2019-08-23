# -*- coding: utf-8 -*-

import re
import h5py
import contextlib
import os
import errno
import re
from time import sleep
import logging
import numpy as np

from . import fs
from . import localfs
from ..utils import instance

logger = logging.getLogger(__name__)


class FileSystemException(fs.FileSystemException):
    """
    Base class for generic HDF5 file system exceptions.
    """
    pass


class LockedError(FileSystemException):
    """
    Device has been locked by someone else.
    """


try:
    unicode
except NameError:
    unicode = str


# Fixed-length dtype's
flen_unicode_dtype = unicode
flen_bytes_dtype = bytes
flen_opaque_dtype = np.void


# Variable-length dtype's
vlen_unicode_dtype = h5py.special_dtype(vlen=unicode)
vlen_bytes_dtype = h5py.special_dtype(vlen=bytes)
vlen_opaque_dtype = h5py.special_dtype(vlen=np.void)


def prepare_h5string(s, raiseExtended=True, useOpaqueDataType=False):
    """
    Convert to Variable-length string (array or scalar).
    Uses UTF-8 encoding when possible.

    Args:
        s: string or array of strings
           string types: unicode, bytes, fixed-length numpy
        raiseExtended(bool): raise decoding error
        useOpaqueDataType(bool): save as opaque instead of
                                 bytes on decoding error
    Returns:
        np.ndarray: dtype is vlen_bytes_dtype, vlen_unicode_dtype
                    or flen_opaque_dtype
    Raises:
        UnicodeDecodeError: extended ASCII encoding
    """
    try:
        # Attempt decoding (not done with vlen_unicode_dtype)
        np.array(s, dtype=flen_unicode_dtype)
    except UnicodeDecodeError:
        # Reason: at least one byte-string with non-ASCII character
        #         (e.g. Latin-1 encoding)
        if raiseExtended:
            raise
        if useOpaqueDataType:
            try:
                np.array(s, dtype=flen_bytes_dtype)
            except UnicodeEncodeError:
                # Reason: at least one UTF-8 string with non-ASCII character
                # Here because: mix of unicode and non-ASCII byte strings
                s = list(map(instance.asbytes, s))
            # vlen_opaque_dtype: currenly not supported by h5py
            return np.array(s, dtype=flen_opaque_dtype)
        else:
            return np.array(s, dtype=vlen_bytes_dtype)
    else:
        return np.array(s, dtype=vlen_unicode_dtype)


def prepare_h5data(data, **kwargs):
    if instance.isstring(data):
        return prepare_h5string(data, **kwargs)
    elif instance.isstringarray(data):
        return prepare_h5string(data, **kwargs)
    else:
        return data


def is_link(fileobj, path):
    """Check whether node is h5py.SoftLink or h5py.ExternalLink

    Args:
        fileobj(File)
        path(str)
    Returns:
        bool
    """
    try:
        lnk = fileobj.get(path, default=None, getlink=True)
    except (KeyError, RuntimeError):
        return False
    else:
        return isinstance(lnk, (h5py.SoftLink, h5py.ExternalLink))


def dereference_link(fileobj, path):
    """Dereference h5py.SoftLink or h5py.ExternalLink
    Args:
        fileobj(File)
        path(str)
    Returns:
        device(str)
        path(str)
    """
    try:
        lnk = fileobj.get(path, default=None, getlink=True)
    except (KeyError, RuntimeError):
        return None, None
    else:
        if isinstance(lnk, h5py.SoftLink):
            #device = fileobj[path].file.filename
            device = fileobj.filename
            dest = lnk.path  # this can be relative to "path"
        elif isinstance(lnk, h5py.ExternalLink):
            device = lnk.filename
            dest = lnk.path
        else:
            device, dest = None, None
        return device, dest


def is_reference(fileobj, path):
    """Check whether node is h5py.Reference

    Args:
        fileobj(File)
        path(str)
    Returns:
        bool
    """
    try:
        lnk = fileobj.get(path, default=None)
    except (KeyError, RuntimeError):
        return False
    else:
        try:
            return h5py.check_dtype(ref=lnk.dtype) == h5py.Reference and lnk.ndim == 0
        except AttributeError:
            return False


def dereference(fileobj, path):
    """Dereference an h5py.Reference

    Args:
        fileobj(File)
        path(str)
    Returns:
        device(str)
        path(str)
    """
    if is_reference(fileobj, path):
        ref = fileobj[path][()]
        dest = fileobj[ref]
        return dest.file.filename, dest.name
    else:
        return None, None


class Enum(dict):
    def __getattr__(self, key):
        return self[key]


h5errnames = ['SIGERR', 'TRUNCERR', 'EEXIST', 'OPENERR']
h5errno = Enum({k: -i for i, k in enumerate(h5errnames, 1)})
h5errcodes = {v: k for k, v in h5errno.items()}


def h5py_errno(err):
    if not errno.errorcode.get(err.errno, None):
        errmsg = str(err)
        m = re.search('errno = ([0-9]+)', errmsg)
        if m:
            return int(m.groups()[0])
        if 'file signature not found' in errmsg:
            return h5errno.SIGERR
        elif 'unable to truncate a file which is already open' in errmsg:
            return h5errno.TRUNCERR
        elif 'file exists' in errmsg:
            return h5errno.EEXIST
        elif 'file is already open for read-only' in errmsg:
            return h5errno.OPENERR
    return err.errno


class h5FileIO(object):
    """Context manager to open HDF5 file with retries
    """

    def __init__(self, path, retries=25, backoff_factor=0.4, **openparams):
        self._handle = None
        self.path = path
        self.retries = retries
        self.backoff_factor = backoff_factor
        self.openparams = openparams

    def __enter__(self):
        self._handle = self._open()
        return self._handle

    def __exit__(self, *args):
        if self._handle:
            self._handle.close()
            self._handle = None

    def _open(self):
        # Open file which is locked by another process:
        #   r   : OK (lock=r), EAGAIN(lock=r+/w/x/a)
        #   r+/w: EAGAIN
        #   x/a : EEXIST
        # Open file which is locked by the same process (possibly another thread):
        #   r : OK
        #   r+: OK(lock=r+/w/x/a), OPENERR(lock=r)
        #   w : TRUNCERR
        #   x : EEXIST
        #   a : EEXIST(lock=r), OK(lock=r+/w/x/a)
        if 'HDF5_USE_FILE_LOCKING' not in os.environ:
            os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
        path = self.path
        for _ in range(self.retries):
            try:
                return h5py.File(path, **self.openparams)
            except IOError as err:
                code = h5py_errno(err)
                if code == errno.ENOENT:
                    raise fs.Missing(path)
                elif code == errno.EISDIR:
                    raise fs.NotAFile(path)
                elif code in [errno.EAGAIN, h5errno.OPENERR, h5errno.TRUNCERR]:
                    pass  # try again
                elif code in [errno.EEXIST, h5errno.EEXIST]:
                    if self.openparams['mode'] in ['x', 'w-']:
                        raise
                    # try again
                else:
                    raise
            sleep(self.backoff_factor)
        raise LockedError(path)


class h5Device(localfs.Path):
    """Proxy to HDF5 file
    """

    def __init__(self, path, mode='a', **kwargs):
        """
        Args:
            path(str):
            mode(Optional(str)):
                r    : open (error when not existing), read
                r+   : open (error when not existing), read/write
                w    : create (truncate existing), read/write
                x/w- : create (error when existing), read/write
                a    : create (append when existing), read/write
        """
        if mode not in ('r', 'r+', 'w', 'x', 'w-', 'a'):
            raise ValueError('Invalid mode {}'.format(repr(mode)))
        self.openparams = kwargs
        self.openparams['mode'] = mode
        super(h5Device, self).__init__(path, **kwargs)

    @contextlib.contextmanager
    def _fopen(self, **openparams):
        #msg = '{} ({})'.format(self.path,openparams)
        #logger.debug('Open '+msg)
        # try:
        with h5FileIO(self.path, **openparams) as f:
            yield f
        # finally:
        #    logger.debug('Close '+msg)

    def _openparams_defaults(self, openparams):
        super(h5Device, self)._openparams_defaults(openparams)
        defaultmode = self.openparams['mode']
        mode = openparams['mode']

        if defaultmode == 'r':
            # read-only
            mode = defaultmode
        elif defaultmode == 'r+' and mode not in ['r', 'r+']:
            # deny new device
            mode = defaultmode
        elif defaultmode in ['x', 'w-'] and mode not in ['r', 'x', 'w-']:
            # allow new device (error when existing)
            mode = defaultmode
        elif defaultmode == 'w' and mode not in ['r', 'w', 'a']:
            # allow new device (overwrite when existing)
            mode = defaultmode
        elif defaultmode == 'a' and mode not in ['r', 'a']:
            # allow new device (append when existing)
            mode = defaultmode

        openparams['mode'] = mode

    @fs.onclose
    def remove(self, **kwargs):
        super(h5Device, self).remove(**kwargs)

    @property
    def mode(self):
        return self.openparams['mode']


class Path(fs.Path):
    """Proxy to HDF5 path
    """

    def __init__(self, path, h5file=None, **kwargs):
        h5file, path = self._split_path(str(path), device=h5file)
        if not path.startswith(self.sep):
            path = self.sep+path
        self.path = path
        if not isinstance(h5file, h5Device):
            if not h5file:
                raise ValueError(
                    'Specify HDF5 file as Path("h5file:/...") or Path("/...",h5file=...)')
            h5file = h5Device(h5file, **kwargs)
        self._h5file = h5file
        super(Path, self).__init__(**kwargs)

    def factory(self, path):
        cls = self.factorycls
        if isinstance(path, cls):
            return path
        else:
            device, path = self._split_path(path, device=self.device)
            return cls(path, h5file=device, **self.factory_kwargs)

    @property
    def openparams(self):
        return self.device.openparams

    @openparams.setter
    def openparams(self, value):
        self.device.openparams = value

    def _openparams_defaults(self, createparams):
        super(Path, self)._openparams_defaults(createparams)
        defaultmode = self.openparams['mode']
        createparams['nodemode'] = createparams['mode']

        # if createparams['nodemode']=='r':
        #    import traceback
        #    print ''
        #    traceback.print_stack()

        # Device access mode (will be further restricted by h5Device)
        if defaultmode == 'w':
            # Do not truncate device when opening files
            createparams['mode'] = 'a'

    @contextlib.contextmanager
    def _fopen(self, **createparams):
        openparams = {k: createparams.pop(k)
                      for k in self.openparams
                      if k in createparams}
        nodemode = createparams.pop('nodemode')

        with self.h5open(**openparams) as f:
            node = f.get(self.path, default=None)
            if node:
                if nodemode == 'w-' or nodemode == 'x':
                    raise fs.AlreadyExists(self.location)
                if nodemode == 'w':
                    if createparams:
                        if isinstance(node, h5py.Group):
                            raise fs.NotAFile(self.location)
                    else:
                        if isinstance(node, h5py.Dataset):
                            raise fs.NotADirectory(self.location)
                    try:
                        del f[self.path]
                    except KeyError:
                        pass
                    node = None
                # r,r+,a: nothing to do
            if node:
                if createparams:
                    raise fs.AlreadyExists(self.location)
            else:
                if nodemode.startswith('r'):
                    raise fs.Missing(self.location)
                parent = self.parent
                node = f.get(parent.path, default=None)
                if not node:
                    raise fs.MissingParentDirectory(self.location)
                node = node.create_dataset(self.name, **createparams)
            yield node

    @contextlib.contextmanager
    def h5open(self, **openparams):
        with self.device.open(**openparams) as f:
            yield f

    def isopen(self):
        return self.device.isopen()

    def h5remove(self):
        self.device.remove()

    @property
    def exists(self):
        # For links: destination exists
        try:
            with self.h5open(mode='r') as f:
                if self.islink:
                    dest = self.linkdest()
                    if dest is None:
                        return False
                    else:
                        return dest.exists
                else:
                    return self.path in f
        except (IOError, fs.FileSystemException):
            return False

    @property
    def lexists(self):
        # For links: link exists
        try:
            with self.h5open(mode='r') as f:
                return self.path in f
        except IOError:
            return False

    @property
    def isdir(self):
        try:
            with self.h5open(mode='r') as f:
                node = f.get(self.path, default=None)
                return isinstance(node, h5py.Group)
        except IOError:
            return False

    @property
    def isfile(self):
        try:
            with self.h5open(mode='r') as f:
                node = f.get(self.path, default=None)
                return isinstance(node, h5py.Dataset)
        except IOError:
            return False

    @property
    def device(self):
        return self._h5file

    @property
    def devsep(self):
        return '::'

    def devsplit(self, path):
        n = len(self.devsep)
        pattern = '\.[^{}]+{}'.format(re.escape(self.sep),
                                      re.escape(self.devsep))
        a = 0
        lst = []
        for m in re.finditer(pattern, path):
            b = m.end()
            add = path[a:b-n]
            if add:
                lst.append(add)
            a = b
        add = path[a:]
        lst.append(add)
        return lst
        # return path.split(self.devsep)

    def listdir(self, recursive=False, depth=0):
        with self.h5open() as f:
            root = f.get(self.path, default=None)
            if isinstance(root, h5py.Group):
                for k in root.keys():
                    yield self[k]

    @property
    def sep(self):
        return '/'

    def mkdir(self, recursive=True, force=True):
        with self.h5open() as f:
            if self._mkdir_prepare(recursive=recursive, force=force):
                return self
            parent = self.parent
            if recursive:
                node = parent
                lst = []
                while not node.exists:
                    lst.append(node)
                    node = node.parent
                for node in reversed(lst):
                    f.create_group(node.path)
            else:
                if not parent.exists:
                    raise fs.MissingParentDirectory(parent.location)
            f.create_group(self.path)
            return self

    def move(self, dest, force=False, rename=False):
        with self.h5open() as f:
            dest = self._copy_move_prepare(dest, force=force, rename=rename)
            if dest.exists and force:
                dest.remove(recursive=True)
            if self.device == dest.device:
                try:
                    # h5py does not have a `rename` so use
                    # the next best thing
                    f.move(self.path, dest.path)
                except ValueError:
                    dest = self._move_copyrenamedel(dest)
            else:
                dest = self._move_copyrenamedel(dest)
            dest._copymove_relink(self, dest)
            return self.factory(dest)

    def _copymove_relink(self, source, dest):
        # Softlink in HDF5 are absolute, so relink all links below self
        # that point to a destination below self (including self)
        for path in self.listdir():
            path._copymove_relink(source, dest)
        if self.islink:
            linkdest = self.linkdest()
            if source.common(linkdest) == source:
                self.remove()
                self.link(dest[source.relpath(linkdest.path)])

    def copy(self, dest, force=False, follow=False, dereference=False):
        with self.h5open() as fsource:
            dest = self._copy_move_prepare(dest, force=force)
            if dest.exists and force:
                dest.remove(recursive=True)
            if self.islink and not follow:
                # just copy the link itself
                dest.link(self.linkdest())
            else:
                with dest.h5open() as fdest:
                    fsource.copy(self.path, fdest[dest.parent.path], name=dest.name,
                                 expand_soft=dereference, expand_external=dereference)
                dest._copymove_relink(self, dest)
            return self.factory(dest)

    def remove(self, recursive=False):
        with self.h5open() as f:
            if self.islink:
                if recursive:
                    node = self.linkdest()
                    if node:
                        node.remove(recursive=True)
                del f[self.path]
            elif self.exists:
                if self.isdir:
                    if not recursive:
                        if list(self.listdir()):
                            raise fs.DirectoryIsNotEmpty(self.location)
                if self.path == self.sep:
                    self.h5remove()
                else:
                    del f[self.path]

    def stats(self, follow=True):
        if follow == False:
            raise ValueError('Hdf5 links do not have attributes themselves')
        with self.open(mode='r') as node:
            return dict(node.attrs)

    def get_stat(self, key, default=None, follow=True):
        if follow == False:
            raise ValueError('Hdf5 links do not have attributes themselves')
        try:
            with self.open(mode='r') as node:
                return node.attrs.get(instance.asunicode(key), default=default)
        except fs.Missing:
            return default

    def pop_stat(self, key, default=None, follow=True):
        if follow == False:
            raise ValueError('Hdf5 links do not have attributes themselves')
        try:
            with self.open(mode='r') as node:
                return node.attrs.pop(instance.asunicode(key), default=default)
        except fs.Missing:
            return default

    def update_stats(self, prepare_kwargs=None, **stats):
        with self.open() as node:
            if prepare_kwargs is None:
                prepare_kwargs = {}
            stats = {instance.asunicode(k): prepare_h5data(v, **prepare_kwargs)
                     for k, v in stats.items()}
            node.attrs.update(stats)

    @contextlib.contextmanager
    def openstats(self, follow=True, **openparams):
        if follow == False:
            raise ValueError('Hdf5 links do not have attributes themselves')
        with self.open(**openparams) as node:
            yield node.attrs

    def link(self, dest, soft=True):
        with self.h5open() as f:
            base = self.parent
            lnkname = self.name

            destpath = self.sibling(dest)
            if destpath.device == self.device:
                if soft:
                    dest = base.relpath(self._getpath(dest))

                    # Remove '..' because it is not supported:
                    # https://support.hdfgroup.org/HDF5/doc/UG/HDF5_Users_Guide-Responsive%20HTML5/index.html#t=HDF5_Users_Guide%2FGroups%2FHDF5_Groups.htm%23IX_soft_links
                    if '..' in dest:
                        # This does not work either:
                        #lst = dest.split(self.sep)
                        #i = len(lst) - lst[::-1].index('..')
                        #base = base[self.sep.join(lst[:i])]
                        #dest = base.relpath(destpath.path)
                        #lnkname = base.relpath(self.path)

                        # Use absolute path instead:
                        dest = destpath.path
                    dest = h5py.SoftLink(dest)
                else:
                    dest = f[destpath.path]
            else:
                dest = h5py.ExternalLink(str(destpath.device), destpath.path)
            f[base.path][lnkname] = dest
            return self.factory(self)

    @property
    def islink(self):
        with self.h5open(mode='r') as f:
            if is_link(f, self.path):
                return True
            else:
                return is_reference(f, self.path)

    def linkdest(self, follow=False):
        if not self.root.exists:
            return None
        with self.h5open(mode='r') as f:
            device, path = dereference_link(f, self.path)
            if path is None:
                device, path = dereference(f, self.path)
            if path and device == self.device:
                dest = self.factory(self.abspath(path))
            elif path and device:
                dest = self.factory(device+self.devsep+path)
            else:
                dest = None
            if follow:
                dest = self._link_follow(dest)
            return dest

    def _contentinfo(self):
        contentinfo = ''
        if self.isfile:
            with self.open(mode='r') as node:
                dtype = node.dtype
                specialdtype = h5py.check_dtype(vlen=node.dtype)
                if node.ndim == 0:
                    if specialdtype == str or specialdtype == unicode:
                        contentinfo = ' = {}'.format(node[()])
                    else:
                        if specialdtype is not None:
                            dtype = specialdtype
                        contentinfo = ' = {} ({})'.format(node[()], dtype)
                else:
                    if specialdtype is not None:
                        dtype = specialdtype
                    shape = 'x'.join(list(map(str, node.shape)))
                    contentinfo = ' = {} ({})'.format(dtype, shape)
        return contentinfo

    def read(self):
        with self.open(mode='r') as node:
            return node[()]

    def write(self, **kwargs):
        return self.mkfile(**kwargs)

    def mkfile(self, data=None, attributes=None, prepare_kwargs=None, **params):
        if prepare_kwargs is None:
            prepare_kwargs = {}
        if data is not None:
            params['data'] = prepare_h5data(data, **prepare_kwargs)
        with self.open(**params) as dset:
            if attributes:
                for k, v in attributes.items():
                    dset.attrs[k] = prepare_h5data(v, **prepare_kwargs)
        return self

    @property
    def dtype(self):
        with self.open(mode='r') as node:
            try:
                return node.dtype
            except AttributeError:
                return None

    @property
    def shape(self):
        with self.open(mode='r') as node:
            try:
                return node.shape
            except AttributeError:
                return None

    @property
    def ndim(self):
        with self.open(mode='r') as node:
            try:
                return node.ndim
            except AttributeError:
                return None

    @property
    def size(self):
        with self.open(mode='r') as node:
            try:
                return node.size
            except AttributeError:
                return None

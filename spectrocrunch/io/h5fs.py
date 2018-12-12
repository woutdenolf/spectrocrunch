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

import re
import h5py
import contextlib
import os
import errno
import re
from time import sleep

from . import fs
from . import localfs

class FileSystemException(fs.FileSystemException):
    """
    Base class for generic HDF5 file system exceptions.
    """
    pass

class LockedError(FileSystemException):
    """
    Device has been locked by someone else.
    """

def h5py_errno(err):
    if errno.errorcode.get(err.errno,None):
        m = re.search('errno = ([0-9]+)',err.message)
        if m:
            return int(m.groups()[0])
    else:
        return err.errno

class h5FileIO(object):
    
    def __init__(self,path,retries=25,backoff_factor=0.4,**openparams):
        self._handle = None
        self.path = path
        self.retries = retries
        self.backoff_factor = backoff_factor
        self.openparams = openparams
    
    def _open(self):
        path = self.path
        for _ in range(self.retries):
            try:
                return h5py.File(path,**self.openparams)
            except IOError as err:
                code = h5py_errno(err)
                if code == errno.ENOENT:
                    raise fs.Missing(path)
                elif code == errno.EISDIR:
                    raise fs.NotAFile(path)
                elif code == errno.EAGAIN:
                    pass
                elif code == errno.EEXIST:
                    if self.openparams['mode'] in ['x','w-']:
                        raise
                else:
                    raise
            sleep(self.backoff_factor)
        raise LockedError(path)

    def __enter__(self):
        self._handle = self._open()
        return self._handle

    def __exit__(self,*args):
        if self._handle:
            self._handle.close()
            self._handle = None

class h5Device(localfs.Path):
    """Proxy to HDF5 file
    """
    
    def __init__(self,path,mode='a',**kwargs):
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
        if mode not in ('r','r+','w','x','w-','a'):
            raise ValueError('Invalid mode {}'.format(repr(mode)))
        self.openparams = kwargs
        self.openparams['mode'] = mode
        super(h5Device,self).__init__(path,**kwargs)
    
    @contextlib.contextmanager
    def _fopen(self,**openparams):
        with h5FileIO(self.path,**openparams) as f:
            yield f
            
    def _openparams_defaults(self,openparams):
        super(h5Device,self)._openparams_defaults(openparams)
        defaultmode = self.openparams['mode']
        mode = openparams['mode']
        if defaultmode=='r':
            # read-only
            mode = defaultmode
        elif defaultmode=='r+' and mode not in ['r','r+']:
            # deny new device
            mode = defaultmode
        elif defaultmode in ['x','w-'] and mode not in ['r','x','w-']:
            # allow new device (error when existing)
            mode = defaultmode
        elif defaultmode=='w' and mode not in ['r','w','a']:
            # allow new device (overwrite when existing)
            mode = defaultmode
        elif defaultmode=='a' and mode not in ['r','a']:
            # allow new device (append when existing)
            mode = defaultmode
        openparams['mode'] = mode
        
    @fs.onclose
    def remove(self,**kwargs):
        super(h5Device,self).remove(**kwargs)
    rm = remove


class Path(fs.Path):
    """Proxy to HDF5 path
    """
    
    def __init__(self,path,h5file=None,**kwargs):
        h5file,path = self._split_path(str(path),device=h5file)
        if not path.startswith(self.sep):
            path = self.sep+path
        self.path = path
        if not isinstance(h5file,h5Device):
            if not h5file:
                raise ValueError('Specify HDF5 file as Path("h5file:/...") or Path("/...",h5file=...)')
            h5file = h5Device(h5file,**kwargs)
        self._h5file = h5file
        super(Path,self).__init__(**kwargs)
        
    @property
    def openparams(self):
        return self.device.openparams

    @openparams.setter
    def openparams(self, value):
        self.device.openparams = value
        
    def _openparams_defaults(self,createparams):
        super(Path,self)._openparams_defaults(createparams)
        defaultmode = self.openparams['mode']
        mode = createparams['mode']
        if defaultmode=='w':
            # Truncate files but not device
            createparams['mode'] = 'a'
        else:
            createparams['mode'] = defaultmode
        createparams['nodemode'] = mode
    
    @contextlib.contextmanager
    def _fopen(self,**createparams):
        openparams = {k:createparams.pop(k)
                      for k in self.openparams 
                      if k in createparams}
        nodemode = createparams.pop('nodemode')

        with self.h5open(**openparams) as f:
            node = f.get(self.path,default=None)
            if node:
                if nodemode=='w-' or nodemode=='x':
                    raise fs.AlreadyExists(self.location)
                if nodemode=='w':
                    if createparams:
                        if isinstance(node,h5py.Group):
                            raise fs.NotAFile(self.location)
                    else:
                        if isinstance(node,h5py.Dataset):
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
                node = f.get(parent.path,default=None)
                if not node:
                    raise fs.MissingParentDirectory(self.location)
                node = node.create_dataset(self.name,**createparams)
            yield node

    @contextlib.contextmanager
    def h5open(self,**openparams):
        with self.device.open(**openparams) as f:
            yield f
    
    def h5remove(self):
        self.device.remove()

    def factory(self,path):
        if isinstance(path,self.__class__):
            return path
        else:
            device,path = self._split_path(path,device=self.device)
            return self.__class__(path,h5file=device,**self.factory_kwargs)

    @property
    def exists(self):
        # For links: destination exists
        try:
            with self.h5open(mode='r') as f:
                if self.islink:
                    return self.linkdest().exists
                else:
                    return self.path in f
        except IOError:
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
                node = f.get(self.path,default=None)
                return isinstance(node,h5py.Group)
        except IOError:
            return False
            
    @property
    def isfile(self):
        try:
            with self.h5open(mode='r') as f:
                node = f.get(self.path,default=None)
                return isinstance(node,h5py.Dataset)
        except IOError:
            return False
    
    @property
    def device(self):
        return self._h5file
    
    @property
    def devsep(self):
        return ':'
        
    def devsplit(self,path):
        n = len(self.devsep)
        pattern = '\.[^{}]+{}'.format(re.escape(self.sep),re.escape(self.devsep))
        a = 0
        lst = []
        for m in re.finditer(pattern,path):
            b = m.end()
            add = path[a:b-n]
            if add:
                lst.append(add)
            a = b
        add = path[a:]
        lst.append(add)
        return lst
        #return path.split(self.devsep)
             
    def listdir(self,recursive=False,depth=0):
        with self.h5open(mode='r') as f:
            root = f.get(self.path,default=None)
            if isinstance(root,h5py.Group):
                for k in root.keys():
                    yield self.factory(self.join(self.path,k))
                
    @property
    def sep(self):
        return '/'

    def mkdir(self, recursive=True, force=True):
        with self.h5open() as f:
            if self._mkdir_prepare(recursive=recursive,force=force):
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
            if self.device==dest.device:
                try:
                    f.move(self.path,dest.path)
                except ValueError:
                    dest = self._move_copydel(dest)
            else:
                dest = self._move_copydel(dest)
            dest._move_relink(self,dest)
            return self.factory(dest)
    
    mv = move
    
    def _move_relink(self,source,dest):
        # Softlink in HDF5 are absolute, so relink all links below self
        # that point to a destination below self (including self)
        for path in self.listdir():
            path._move_relink(source,dest)
        if self.islink:
            linkdest = self.linkdest()
            if source.common(linkdest)==source:
                self.remove()
                self.link(dest[source.relpath(linkdest.path)])

    def copy(self, dest, force=False, follow=False, dereference=False):
        with self.h5open() as fsource:
            dest = self._copy_move_prepare(dest, force=force)
            if dest.exists and force:
                dest.remove(recursive=True)
            if self.islink and not follow:
                # just copy the link
                dest.link(self.linkdest)
            else:
                # TODO: missing option: expand softlinks when pointing outside
                #       the tree being copied.
                with dest.h5open() as fdest:
                    fsource.copy(self.path,fdest[dest.parent.path],name=dest.name,
                                 expand_soft=dereference,expand_external=dereference)
            return self.factory(dest)
    
    cp = copy
    
    def remove(self,recursive=False):
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
                if self.path==self.sep:
                    self.h5remove()
                else:
                    del f[self.path]
    
    rm = remove
    
    def stats(self,follow=True):
        if follow==False:
            raise ValueError('Hdf5 links do not have attributes themselves')
        with self.open(mode='r') as node:
            return dict(node.attrs)

    def get_stat(self,key,default=None,follow=True):
        if follow==False:
            raise ValueError('Hdf5 links do not have attributes themselves')
        try:
            with self.open(mode='r') as node:
                return node.attrs.get(key,default=default)
        except fs.Missing:
            return default
            
    def update_stats(self,**stats):
        with self.open() as node:
            node.attrs.update(stats)
    
    @contextlib.contextmanager
    def openstats(self,follow=True,**openparams):
        if follow==False:
            raise ValueError('Hdf5 links do not have attributes themselves')
        with self.open(**openparams) as node:
            yield node.attrs
            
    def link(self,dest,soft=True):
        with self.h5open() as f:
            base = self.parent
            lnkname = self.name
            
            destpath = self.sibling(dest)
            if destpath.device==self.device:
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
            try:
                lnk = f.get(self.path,default=None,getlink=True)
            except KeyError:
                lnk = None
            return isinstance(lnk,(h5py.SoftLink,h5py.ExternalLink))

    def linkdest(self,follow=False):
        if not self.root.exists:
            return None
        with self.h5open(mode='r') as f:
            try:
                lnk = f.get(self.path,default=None,getlink=True)
            except KeyError:
                lnk = None
            else:
                if isinstance(lnk,h5py.SoftLink):
                    # TODO: h5py.SoftLink does not contain information on the file
                    #       which can cause problems when some of the parents are
                    #       external links. Should check for that here.
                    lnk = self.factory(self.abspath(lnk.path))
                elif isinstance(lnk,h5py.ExternalLink):
                    lnk = self.factory("{}{}{}".format(lnk.filename,self.devsep,lnk.path))
                else:
                    lnk = None
            if follow:
                lnk = self._link_follow(lnk)
            return lnk
    
    def _contentinfo(self):
        contentinfo = ''
        if self.isfile:
            with self.open(mode='r') as node:
                dtype = node.dtype
                specialdtype = h5py.check_dtype(vlen=node.dtype)
                if node.ndim==0:
                    if specialdtype==str or specialdtype==unicode:
                        contentinfo = ' = {}'.format(node[()])
                    else:
                        if specialdtype is not None:
                            dtype = specialdtype
                        contentinfo = ' = {} ({})'.format(node[()],dtype)
                else:
                    if specialdtype is not None:
                        dtype = specialdtype
                    shape = 'x'.join(list(map(str,node.shape)))
                    contentinfo = ' = {} ({})'.format(dtype,shape)
        return contentinfo

    def read(self):
        with self.open(mode='r') as node:
            return node[()]
    
    def write(self,**kwargs):
        return self.mkfile(**kwargs)
    

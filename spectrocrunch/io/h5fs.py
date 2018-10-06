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

from . import fs

class h5File(fs.File):

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
            raise ValueError('Invalid mode \'{}\''.format(mode))
        self.openparams = kwargs
        self.openparams['mode'] = mode
        self.path = path
        super(h5File,self).__init__(**kwargs)
    
    @contextlib.contextmanager
    def _fopen(self,**openparams):
        with h5py.File(self.path,**openparams) as f:
            yield f

    @property
    def openparams(self):
        return self._openparams
    
    @openparams.setter
    def openparams(self,value):
        self._openparams = value
        
        
class Path(fs.Path):

    def __init__(self,path,h5file=None,**kwargs):
        h5file,self.path = self._split_path(str(path),device=h5file)
        if not isinstance(h5file,h5File):
            if not h5file:
                raise ValueError('Specify HDF5 file as Path("h5file:/...") or Path("/...",h5file=...)')
            h5file = h5File(h5file,**kwargs)
        self._h5file = h5file
        super(Path,self).__init__(**kwargs)
    
    @property
    def openparams(self):
        return self._h5file.openparams

    @contextlib.contextmanager
    def _fopen(self,**createparams):
        keys = (set(self.openparams.keys()) & set(createparams.keys()))
        openparams = {k:createparams.pop(k) for k in keys}
        mode = openparams.get('mode',self.openparams.get('mode','a'))

        with self.h5open(**openparams) as f:
            node = f.get(self.path,default=None)
            if node:
                if mode=='w-' or mode=='x':
                    raise fs.AlreadyExists(self.location)
            
                if mode=='w':
                    if createparams:
                        if isinstance(node,h5py.Group):
                            raise fs.NotAFile(self.location)
                    else:
                        if isinstance(node,h5py.Dataset):
                            raise fs.NotADirectory(self.location)
                    del f[self.path]
                    node = None
                # r,r+,a: nothing to do
                
            if not node:
                if mode.startswith('r'):
                    raise fs.Missing(self.location)
                
                parent = self.parent
                node = f.get(parent,default=None)
                if not node:
                    raise fs.MissingParentDirectory(self.location)
                node = node.create_dataset(self.name,**createparams)
            yield node

    @contextlib.contextmanager
    def h5open(self,**openparams):
        # Do not truncate (w) or exclusive create (x/w-)
        mode = openparams.get('mode',self.openparams.get('mode','a'))
        if mode.startswith('w') or mode=='x':
            mode = 'a'
        openparams['mode'] = mode
        with self._h5file.open(**openparams) as f:
            yield f
    
    def __eq__(self,other):
        return self.location==str(getattr(other,'location',other))

    def __ne__(self,other):
        return self.location!=str(getattr(other,'location',other))
        
    def factory(self,path):
        if isinstance(path,self.__class__):
            return path
        else:
            device,path = self._split_path(path,device=self._h5file)
            return self.__class__(path,h5file=device)

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
        return self._h5file.path
    
    def _split_path(self,path,device=None):
        """
        Args:
            path(str):
            device(Optional(str|h5File)): when device description in path is missing
            
        Returns:
            device(str):
            path(str):
        """
        device,path = super(Path,self)._split_path(path,device=device)
        if not path:
            path = self.sep
        return device,path
             
    def norm(self,path):
        path = self._split_path(path)[1]
        if not path.startswith(self.sep):
            path = self.sep+path
        return super(Path,self).norm(path)
        
    def listdir(self,recursive=False,depth=0):
        with self.h5open(mode='r') as f:
            root = f.get(self.path,default=None)
            if isinstance(root,h5py.Group):
                for k in root.keys():
                    yield self.factory(self.join(self.path,k))
            elif root:
                yield root.name
                
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
            
    def move(self, dest, force=True):
        with self.h5open() as f:
            dest = self._copy_move_prepare(dest, force=force)
            f.move(self.path,dest.path)
            return dest
        
    def copy(self, dest, force=True, follow=False):
        with self.h5open() as fsource:
            dest = self._copy_move_prepare(dest, force=force)
            with dest.h5open() as fdest:
                fsource.copy(self.path,fdest[dest.parent],name=dest.name,
                             expand_soft=follow,expand_external=follow)
                return dest
        
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
                del f[self.path]
                
    def stats(self,follow=True):
        with self.open(mode='r') as node:
            return dict(node.attrs)

    def update_stats(self,**stats):
        with self.open() as node:
            node.attrs.update(stats)
        
    def link(self,dest,soft=True):
        with self.h5open() as f:
            lnkname = self.path
            _dest = self.factory(dest)
            if _dest.device!=self.device:
                f[lnkname] = h5py.ExternalLink(_dest.device, _dest.path)
            else:
                if soft:
                    dest = self.relpath(str(dest))
                    f[lnkname] = h5py.SoftLink(dest)
                else:
                    dest = str(_dest)
                    f[lnkname] = f[dest]
            return self
            
    @property
    def islink(self):
        with self.h5open(mode='r') as f:
            try:
                lnk = f.get(self.path,default=None,getlink=True)
            except KeyError:
                lnk = None
            return isinstance(lnk,(h5py.SoftLink,h5py.ExternalLink))

    def linkdest(self,follow=False):
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
        

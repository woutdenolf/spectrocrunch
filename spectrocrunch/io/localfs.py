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
import errno
import shutil
import contextlib
from datetime import datetime

from . import fs

class Path(fs.Path):

    def __init__(self,path,mode='a',**kwargs):
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
        if mode not in ('r','r+','w','w+','x','x+','a','a+'):
            raise ValueError('Invalid mode \'{}\''.format(mode))
        self.openparams = kwargs
        self.openparams['mode'] = mode
        
        self.path = path
        super(Path,self).__init__(**kwargs)
    
    @property
    def openparams(self):
        return self._openparams
    
    @openparams.setter
    def openparams(self,value):
        self._openparams = value
    
    @contextlib.contextmanager
    def _fopen(self,**openparams):
        try:
            with open(self.path,**openparams) as f:
                yield f
        except IOError as err:
            if err.errno == errno.ENOENT:
                raise fs.Missing(self.location)
            elif err.errno == errno.EISDIR:
                raise fs.NotAFile(self.location)
            else:
                raise 
 
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
        for k in os.listdir(self.path):
            yield self.factory(self.join(self.path,k))

    @property
    def sep(self):
        return os.sep

    def mkdir(self,recursive=True,force=True):
        if self._mkdir_prepare(recursive=recursive,force=force):
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
        
    def move(self, dest, force=True):
        dest = self._copy_move_prepare(dest, force=force)
        try:
            os.rename(self.path, dest.path)
        except OSError as err:
            if err.errno == errno.EXDEV:
                dest = self._move_copydel(dest)
            else:
                raise
                
        return dest
                
    def copy(self, dest, force=True, follow=False):
        dest = self._copy_move_prepare(dest, force=force)
        if self.isdir:
            shutil.copytree(self.path, dest.path, symlinks=not follow)
        else:
            shutil.copy(self.path, dest.path)
        return dest

    def remove(self,recursive=False):
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
                    if err.errno==errno.ENOTEMPTY:
                        raise fs.DirectoryIsNotEmpty(self.location)
                    else:
                        raise
        elif self.exists:
            os.remove(self.path)
        
    def stats(self,follow=True):
        if follow:
            ret = os.lstat(self.path)
        else:
            ret = os.stat(self.path)
        ret = {k[3:]:getattr(ret,k) for k in dir(ret) if k.startswith('st_')}
        
        for k,v in ret.items():
            if 'time' in k:
                ret[k] = datetime.fromtimestamp(v)
                
        return ret
        
    def update_stats(self,permissions=None,uid=-1,gid=-1):
        if permissions is not None:
            os.chmod(self.path,permissions)
        if uid>0 or gid>0:
            os.chown(self.path,uid,gid)
        
    def link(self,dest,soft=True):
        lnkname = self.path
        if soft:
            dest = self.relpath(str(dest))
            os.symlink(dest, lnkname)
        else:
            dest = str(self.factory(dest))
            os.link(dest, lnkname)
        return self
        
    def linkdest(self,follow=False):
        try:
            lnk = self.factory(self.abspath(os.readlink(self.path)))
            if follow:
                lnk = self._link_follow(lnk)
            return lnk
        except OSError:
            return None
            

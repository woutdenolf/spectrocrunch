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
from abc import ABCMeta,abstractproperty,abstractmethod
from future.utils import with_metaclass
import operator

from . import utils

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
    
    
class File(with_metaclass(ABCMeta,object)):

    def __init__(self,**kwargs):
        self._handle = None
        self._onclose_callbacks = []
        super(File,self).__init__()

    @contextlib.contextmanager
    def open(self,**openparams):
        if self._handle is None:
            kwargs = dict(self.openparams)
            kwargs.update(openparams)
            with self._closectx():
                with self._fopen(**kwargs) as self._handle:
                    try:
                        yield self._handle
                    finally:
                        self._handle = None
        else:
            yield self._handle

    @contextlib.contextmanager
    @abstractmethod
    def _fopen(self,**openparams):
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
                funcs,self._onclose_callbacks = self._onclose_callbacks,[]
                for func,args,kwargs in funcs:
                    func(*args,**kwargs)
    
    def executeonclose(self,func,*args,**kwargs):
        if self._handle:
            self._onclose_callbacks.append((func,args,kwargs))
        
    def isopen(self):
        # only checks local handle (could be opened by others)
        return bool(self._handle)
        
        
class Path(File):

    @property
    def path(self):
        return self._path
    
    @path.setter
    def path(self,value):
        self._path = self.norm(value)
    
    @property
    def location(self):
        if self.device:
            return '{}{}{}'.format(self.device,self.devsep,self.path)
        else:
            return self.path
            
    @property
    def device(self):
        return ''
    
    @property
    def devsep(self):
        return ':'
    
    def _split_path(self,path,device=None):
        """
        Args:
            path(str):
            device(Optional(str)): when device description in path is missing
            
        Returns:
            device(str):
            path(str):
        """
        tmp = str(path).split(self.devsep)
        n = len(tmp)
        if n==2:
            device,path = tmp
        elif n>2:
            device,path = tmp[0],self.devsep.join(tmp[1:])
        else:
            path = tmp[0]
        return device,path

    @classmethod
    def factory(cls,path):
        if isinstance(path,cls):
            return path
        else:
            return cls(path)
            
    def __repr__(self):
        return "'{}'".format(self.path)
        
    def __str__(self):
        return self.path
    
    def __hash__(self):
        return hash(str(self))

    def __bool__(self):
        return bool(str(self))
        
    def __nonzero__(self):
        return self.__bool__()
    
    def _binary_op(self,other,op):
        return op(str(self),str(other))

    def __eq__(self,other):
        return self._binary_op(other,operator.eq)
    
    def __ne__(self,other):
        return self._binary_op(other,operator.ne)
    
    def __le__(self,other):
        return self._binary_op(other,operator.le)
    
    def __ge__(self,other):
        return self._binary_op(other,operator.ge)
        
    def __lt__(self,other):
        return self._binary_op(other,operator.lt)
    
    def __gt__(self,other):
        return self._binary_op(other,operator.gt)
    
    def __len__(self):
        return len(str(self))
        
    def __add__(self,other):
        return self[other]

    def __iadd__(self,other):
        self.path += str(other)

    def split(self,*args):
        return str(self).split(self.sep)
    
    def endswith(self,*args):
        return str(self).endswith(*args)
    
    def startswith(self,*args):
        return str(self).endswith(*args)
    
    def encode(self,*args):
        return str(self).encode(*args)
        
    def __getitem__(self,value):
        return self.factory(self.join(self.path,str(value)))

    def __call__(self,*value):
        return self.factory(self.join(self.path,*map(str,value)))

    def append(self,value):
        self += value

    @property
    def isfile(self):
        return not self.isdir

    @property
    def parent(self):
        return self['..']
            
    def _copy_move_prepare(self, dest, force=True):
        dest = self.factory(dest)
        if dest.exists and not force:
            raise AlreadyExists(dest.location)
            
        d = dest.parent
        if d and not d.exists:
            d.mkdir()
        
        return dest

    def _move_copydel(self,dest):
        """This is called if dest is an external path
           and cannot be move directly. So copy/delete
           instead, which is atomic if rename is atomic.
        """
        tmp = '{}-{}'.format(dest.name, utils.randomstring())
        tmp = dest.parent[tmp]
        self.copy(tmp,force=True,follow=True)
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

    def _sep_in(self,path):
        if os.sep==self.sep:
            return path
        else:
            return path.replace(self.sep,op.sep)

    def _sep_out(self,path):
        if os.sep==self.sep:
            return path
        else:
            return path.replace(os.sep,self.sep)
            
    def join(self,*args):
        args = list(map(self._sep_in,args))
        return self._sep_out(os.path.join(*args))

    def norm(self,path):
        path = self._sep_in(path)
        path = os.path.abspath(os.path.normpath(path))
        return self._sep_out(path)
        
    def abspath(self,path):
        return self.norm(self.join(self.parent.path,path))
        
    def relpath(self,path):
        path = self._sep_in(path)
        path = os.path.relpath(path,self.parent.path)
        return self._sep_out(path)

    @property
    def name(self):
        return os.path.basename(self._sep_in(self.path))
        
    @property
    def root(self):
        root = self
        while root!=root.parent and root.parent:
            root = root.parent
        return root
        
    def tree(self,recursive=True,depth=0,files=True,_level=0):
        """
        Returns:
            2-tuple: (Path|None, list(2-tuple)|None)
        """
        if self.isdir:
            lst = []
            out = (self,lst)

            for path in sorted(self.listdir()):
                if path.isdir:
                    if not recursive and _level==depth:
                        lst.append((path,None))
                    else:
                        lst.append(path.tree(recursive=recursive,depth=depth,
                                             files=files,_level=_level+1))
                elif files:
                    lst.append((path,None))
             
        else:
            if files:
                out = (self,None)
            else:
                out = (None,None)
        
        return out

    def strtree(self,recursive=True,depth=0,files=True,stats=False):
        tree = self.tree(recursive=recursive,depth=depth,files=files)
        if stats:
            lines = self._str_tree_withstats(tree)
        else:
            lines = self._str_tree(tree)
        return '\n'.join(lines)
        
    def ls(self,recursive=True,depth=0,files=True,stats=False):
        tree = self.strtree(recursive=recursive,depth=depth,files=files,stats=stats)
        print(tree)

    @classmethod
    def _str_tree(cls,node,_padding=' '):
        path,lst = node
        out = []
        if path:
            # Node name
            indent = _padding[:-1]
            
            if indent:
                node = path.name
            else:
                node = path.location
            nodename = '{}+-{}'.format(indent,node)
            
            # Add link info
            lnkdest = path.linkdest(follow=False)
            if lnkdest:
                if lnkdest.lexists:
                    sep = '-'
                else:
                    sep = '/'
                nodename = '{} -{}-> '.format(nodename,sep)
                
                if path.device == lnkdest.device:
                    lnkdest = path.relpath(str(lnkdest))
                    lst = None
                else:
                    lnkdest = lnkdest.location
                    _padding += ' '*(len(nodename)-4)
                    
                nodename += lnkdest

            out.append(nodename)
            
            # Add subnodes
            if lst:
                last = len(lst)
                _padding = _padding + ' '
                for i,node in enumerate(lst,1):
                    out.append(_padding+'|')
                    if i==last:
                        out += cls._str_tree(node,_padding=_padding+' ')
                    else:
                        out += cls._str_tree(node,_padding=_padding+'|')
        return out

    @classmethod
    def _str_tree_withstats(cls,node,_padding=' '):
        path,lst = node
        out = []
        if path:
            # Node name
            indent = _padding[:-1]
            name,stats = path.lsinfo()
            out.append('{}+-{}'.format(indent,name))
            for k,v in stats.items():
                out.append('{} | @{} = {}'.format(_padding,k,v))
                
            # Add subnodes
            if lst:
                last = len(lst)
                _padding = _padding + ' '
                for i,node in enumerate(lst,1):
                    out.append(_padding+'|')
                    if i==last:
                        out += cls._str_tree_withstats(node,_padding=_padding+' ')
                    else:
                        out += cls._str_tree_withstats(node,_padding=_padding+'|')
        return out
    
    def lsinfo(self):
        try:
            stats = self.stats()
            name = self.name
        except Missing:
            name = '{} (broken link: {})'.format(self.name,self.linkdestname)
            stats = {}
        return name,stats
        
    def node(self,follow=True):
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

    @abstractmethod
    def move(self, dest, force=True):
        pass

    def rename(self, dest):
        return self.move(dest, force=False)
        
    @abstractmethod
    def copy(self, dest, force=True, follow=False):
        pass

    @abstractmethod
    def stats(self,follow=True):
        pass
    
    @abstractmethod
    def update_stats(self,follow=True,**stats):
        pass
    
    @contextlib.contextmanager
    def openstats(self,follow=True):
        stats = self.stats(follow=follow)
        yield stats
        self.update_stats(**stats)
    
    @abstractmethod
    def remove(self,recursive=False):
        pass
    
    @abstractmethod
    def link(self,dest,soft=True):
        pass
    
    @abstractmethod
    def linkdest(self,follow=False):
        pass
        
    @property
    def linkdestname(self):
        lnkdest = self.linkdest(follow=False)
        if lnkdest:
            if self.device == lnkdest.device:
                lnkdest = self.relpath(str(lnkdest))
            else:
                lnkdest = lnkdest.location
        return lnkdest

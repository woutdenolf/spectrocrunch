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

import collections
import itertools
import contextlib
try:
    import itertools.imap as map
except ImportError:
    pass

from ..common import hashing
from ..common import listtools

class AnnotatedKey(object):

    def __init__(self,key,localcontext=False,changed=False):
        self.original = key
        self.localcontext = localcontext
        self.changed = changed
        
    @property
    def inherited(self):
        return not self.local
    
    def __repr__(self):
        if self.localcontext:
            s = ""
        else:
            s = "I"
        if self.changed:
            s = s+"*"
        if s:
            return "{}({})".format(self.original,s)
        else:
            return str(self.original)
    
class Node(collections.MutableMapping):
    # Original idea:
    # http://code.activestate.com/recipes/577434-nested-contexts-a-chain-of-mapping-objects/

    def __init__(self, parent=None, name=None):
        self._name = name
        
        # Local parameters
        self._local_map = {}
        self.rehash()
        self.changed = True
        
        # Link to other nodes
        self._parent = None
        self.children = []
        self.parent = parent

    @property
    def name(self):
        if self._name:
            return self._name
        else:
            return str(id(self))

    @property
    def parent(self):
        return self._parent
    
    @parent.setter
    def parent(self,node):
        if node is self or node is self.parent:
            return
            
        if node is not None:
            # "node" is downstream from "self:"
            # make "node" a sibling to "self"
            if id(node) in map(id,self.iter_down):
                node._newparent(self.parent)
                
        self._newparent(node)
        
    def _newparent(self,node):
        # Remove self from self.parent:
        if self.parent is not None:
            self.parent.children = [child for child in self.parent.children if child is not self]
        # Change self.parent:
        self._parent = node
        # 
        if node is not None:
            node.children.append(self)
            
    def branch(self, **kwargs):
        """Add a node inheriting from this one
        """
        return self.__class__(parent=self,**kwargs)

    def insert_before(self,node=None,single=False,**kwargs):
        """Insert a node before this one
        """
        newparent = None
        if node is self.parent or node is self:
            return node,newparent
        if node is None:
            node = self.__class__(**kwargs)
        if single:
            newparent = node.remove(single=True)
        self.parent,keep = node,self.parent
        if node.parent is None:
            node.parent = keep
        return node,newparent
        
    def insert_after(self,node=None,single=False,**kwargs):
        """Insert a node after this one
        """
        newparent = None
        if self is node.parent or node is self:
            return node,newparent
        if node is None:
            node = self.__class__(**kwargs)
        if single:
            newparent = node.remove(single=True)
        for child in list(self.children):
            child.parent = node
        node.parent = self
        return node,newparent
    
    def remove(self,single=True):
        newparent = None
        if single:
            # Move children to parent
            if self.parent is None:
                newparent = self.__class__()
                self.parent = newparent
            parent = self.parent
            for child in list(self.children):
                child.parent = parent
            if parent is not self.parent:
                print parent,self.parent
        self.parent = None
        return newparent
        
    @property
    def dtype(self):
        return self._local_map.__class__
        
    def todict(self,full=True,annotated=False):
        if full:
            return self.dtype(self.items(annotated=annotated))
        else:
            return self.dtype(self._map(annotated=annotated))
            
    def fromdict(self,dic,override=True,recursive=False):
        if not override:
            dic = {k:v for k,v in dic.items() if not self.inherited(k)}
    
        if recursive:
            self.update_recursive(dic)
        else:
            self.update(dic)
    
    def update_recursive(self,*args,**kwargs):
        """Update map locally and remove _children override of the keys
        """
        d = self.dtype(*args, **kwargs)
        self.update(d)
        for key in d:
            self.delete_recursive(key,local=False)
              
    @property
    def _root(self):
        return self if self.parent is None else self.parent._root

    def _annotate_key(self,key,localcontext=True):
        return AnnotatedKey(key,localcontext,changed=key in self._newkeys)
        
    def _map(self,annotated=False,localcontext=True):
        # Returns a copy when annotated
        if annotated:
            if self.dtype is dict:
                return {self._annotate_key(k,localcontext=localcontext):v for k,v in self._local_map.items()}
            else:
                return self.dtype(map(lambda item:(self._annotate_key(item[0],localcontext=localcontext),item[1]),self._local_map.items()))
        else:
            return self._local_map

    def _maps(self,annotated=False):
        """Generator of all node maps in ascending order of inheritance (self first)
        """
        yield self._map(annotated=annotated)
        node = self.parent
        while node is not None:
            yield node._map(annotated=annotated,localcontext=False)
            node = node.parent

    @property
    def iter_up(self):
        """Generator of all nodes in ascending order (self first)
        """
        yield self
        node = self.parent
        while node is not None:
            yield node
            node = node.parent

    @property
    def iter_down(self):
        yield self
        for child in self.children:
            for node in child.iter_down:
                yield node
    
    @property
    def root(self):
        for node in self.iter_up:
            pass
        return node

    def __getitem__(self, key):
        for m in self._maps():
            if key in m:
                break
        return m[key]

    @contextlib.contextmanager
    def _changekey(self,key):
        try:
            yield self._map()
        except Exception as e:
            raise e
        else:
            # Key has changed
            self._newkeys.add(key)
            self.rehash()

    def rehash(self,recursive=False):
        self._hash_map_cached = None
        if recursive: # not sure we'll even need it
            for child in self.children:
                child.rehash(recursive=recursive)

    def __setitem__(self, key, value):
        try:
            if hashing.hashequal(self[key],value):
                return
        except KeyError:
            pass
        
        with self._changekey(key) as cmap:
            cmap[key] = value
        
    def delete_local(self,key):
        """Delete a key locally
        """
        try:
            with self._changekey(key) as cmap:
                del cmap[key]
        except KeyError:
            pass
            
    def delete_recursive(self,key,local=True):
        """Delete a key in all _children and optionally locally
        """
        if local:
            self.delete_local(key)
        for child in self.children:
            child.delete_recursive(key)
    
    def __delitem__(self, key):
        """Delete a key locally and in all _children
        """
        self.delete_recursive(key)

    def __len__(self):
        return len(set(itertools.chain.from_iterable(self._maps())))

    def __iter__(self):
        return self.iter(annotated=False)

    def iter(self,annotated=False):
        # Iterate over keys
        it = itertools.chain.from_iterable(self._maps(annotated=annotated))
        if annotated:
            key = lambda key:key.original
        else:
            key = None
        return listtools.unique_everseen(it,key=key)
        
    def keys(self,annotated=False):
        return self.iter(annotated=annotated)
    
    def items(self,annotated=False):
        if annotated:
            return map(lambda key:(key,self[key.original]),self.keys(annotated=True))
        else:
            return super(Node,self).items()
            
    def __contains__(self, key):
        return any(key in m for m in self._maps())

    def inherited(self, key):
        return any(key in m for m in self._maps()[1:])

    def tree(self,level=0,onlychanged=False):
        """Show inheritance tree of this node
        """
        prefix = "  "*level
        if prefix:
            prefix += "> "

        lst = [child.tree(level=level+1,onlychanged=onlychanged) for child in self.children]
        lst.insert(0,"{}{} = {}".format(prefix,self.name,self._repr_local(onlychanged=onlychanged)))
        
        return "\n".join(lst)
        
    def _repr_local(self, onlychanged=False):
        if onlychanged:
            lst = ["{}: {}".format(k,v) for k,v in self._map().items() if k in self._newkeys]
        else:
            lst = ["{}{}: {}".format(k,"(*)" if k in self._newkeys else "",v) for k,v in self._map().items()]
            
        if self.changed:
            return "{"+", ".join(lst)+"}(*)"
        else:
            return "{"+", ".join(lst)+"}"

    def _repr_all(self, onlychanged=False):
        lst = ["{}: {}".format(k,v) for k,v in self.items(annotated=True)]
        if self.changed:
            return "{"+", ".join(lst)+"}(*)"
        else:
            return "{"+", ".join(lst)+"}"

    def _repr_inheritance(self):
        return  " -> ".join([node._repr_local() for node in self.iter_up])

    def __repr__(self):
        return "{} = {}".format(self.name,self._repr_all())

    @property
    def _hash_local(self):
        if not self._hash_map_cached:
            self._hash_map_cached = hashing.calchash(self._map())
        return self._hash_map_cached

    def __hash__(self):
        # Hash of all locals up the tree
        return hash((self._hash_local,hash(self.parent)))
  
    @property
    def changed(self):
        return hash(self)!=self._hash_all_cached

    @changed.setter
    def changed(self,changedstate):
        if changedstate:
            self._hash_all_cached = None
            self._newkeys = set(self._map())
        else:
            self._hash_all_cached = hash(self)
            self._newkeys = set()
        
    def reset(self,recursive=False,changedstate=False):
        self.changed = changedstate
        if recursive:
            for child in self.children:
                child.reset(recursive=recursive,changedstate=changedstate)


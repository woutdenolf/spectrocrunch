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
#from ordered_set import OrderedSet
from ..common import hashing
from ..common import listtools

class Node(collections.MutableMapping):
    # Original idea:
    # http://code.activestate.com/recipes/577434-nested-contexts-a-chain-of-mapping-objects/

    def __init__(self, parent=None, name=None, dtype=dict):
        self._name = name
        
        # Data storage:
        if parent:
            cls = parent.dtype
        else:
            cls = dtype
        self._map = cls()
        self._rehash = True
        self._verify_hash = None
        
        # Link with other nodes:
        self._parent = parent
        self._children = []
        if parent:
            parent._children.append(self)

    @property
    def name(self):
        if self._name:
            return self._name
        else:
            return str(id(self))

    def branch(self, **kwargs):
        """Add a node inheriting from this one
        """
        return self.__class__(parent=self,**kwargs)
        
    def insert_before(self,**kwargs):
        """Insert a node before this one
        """
        parent = self._parent
        
        # Add child to parent
        child = self.__class__(parent=parent,**kwargs)
        
        # Link self with child
        self._parent = child
        child._children.append(self)
        
        # Remove self from parent
        if parent:
            parent._children.remove(self)

        return child
        
    def insert_after(self,**kwargs):
        """Insert a node after this one
        """
        # Add child to self
        child = self.__class__(parent=self,**kwargs)
        
        # Link self._children with child
        print self._children
        for child2 in self._children:
            print child2,child,child2==child
            if child2 != child:
                child2._parent = child
                child._children.append(child2)
                
        # Keep only child in self._children
        self._children = [child]
        
        return child
    
    def tree(self,level=0):
        """Show inheritance tree of this node
        """
        ret = "  "*level
        if ret:
            ret += "> "

        ret = "{}{} = {}\n".format(ret,self.name,self._map)
        for child in self._children:
            ret += child.tree(level=level+1)
            
        return ret
    
    @property
    def dtype(self):
        return self._map.__class__
        
    def todict(self,full=True):
        if full:
            return self.dtype(self.items())
        else:
            return self.dtype(self._map)
    
    def fromdict(self,dic,override=True,recursive=False):
        if not override:
            dic = {k:v for k,v in dic.items() if not self.inherited(k)}
    
        if recursive:
            self.update_recursive(dic)
        else:
            self.update(dic)
    
    def update_recursive(self,*args,**kwargs):
        """Update map locally and remove children override of the keys
        """
        d = self.dtype(*args, **kwargs)
        self.update(d)
        for key in d:
            self.delete_recursive(key,local=False)
              
    @property
    def _root(self):
        return self if self._parent is None else self._parent._root

    @property
    def _maps(self):
        """Generator of all node maps in ascending order of inheritance (self first)
        """
        yield self._map
        node = self._parent
        while node:
            yield node._map
            node = node._parent
    
    @property
    def _nodes_up(self):
        """Generator of all nodes in ascending order of inheritance (self first)
        """
        yield self
        node = self._parent
        while node:
            yield node
            node = node._parent
            
    @property
    def _rmaps(self):
        """Iterator of all node maps in desceding order of inheritance (root first)
        """
        return reversed(list(self._maps))
    
    @property
    def _nodes_rup(self):
        """Iterator of all nodes in desceding order of inheritance (root first)
        """
        return reversed(list(self._nodes_up))
    
    def __getitem__(self, key):
        for m in self._maps:
            if key in m:
                break
        return m[key]

    @property
    def _cmap(self):
        self.rehash(recursive=False)
        return self._map

    def rehash(self,recursive=True):
        self._rehash = True
        if recursive:
            for child in self._children:
                child.rehash(recursive=recursive)

    def __setitem__(self, key, value):
        try:
            if hashing.hashequal(self[key],value):
                return
        except KeyError:
            pass
                 
        self._cmap[key] = value
        
    def delete_local(self,key):
        """Delete a key locally
        """
        try:
            del self._cmap[key]
        except KeyError:
            pass
            
    def delete_recursive(self,key,local=True):
        """Delete a key in all children and optionally locally
        """
        if local:
            self.delete_local(key)
        for child in self._children:
            child.delete_recursive(key)
    
    def __delitem__(self, key):
        """Delete a key locally and in all children
        """
        self.delete_recursive(key)

    def __len__(self):
        return len(set(itertools.chain.from_iterable(self._maps)))

    def __iter__(self):
        return listtools.unique_everseen(itertools.chain.from_iterable(self._maps))
        #return iter(OrderedSet(itertools.chain.from_iterable(self._maps)))
        #return iter(set(itertools.chain.from_iterable(self._maps)))

    def __contains__(self, key):
        return any(key in m for m in self._maps)

    def inherited(self, key):
        return any(key in m for m in self._maps[1:])

    def __repr__(self):
        s = " <- ".join([repr(node._map)+("(*)" if node.changed else "") for node in self._nodes_rup])
        return "{} = {}".format(self.name,s)

    @property
    def _hash_local(self):
        if self._rehash:
            hself = hashing.calchash(self._map)
            self._hself = hself
            self._rehash = False
        else:
            hself = self._hself
        return hself

    def __hash__(self):
        return hash((self._hash_local,hash(self._parent)))

    def __eq__(self,other):
        if isinstance(other,self.__class__):
            return hash(self)==hash(other)
        else:
            return hashing.hashequal(self.todict(),other)
            
    def __neq__(self,other):
        return not self.__eq__(other)
  
    @property
    def changed(self):
        return hash(self)!=self._verify_hash

    @changed.setter
    def changed(self,changedstate):
        hself = hash(self)
        if changedstate:
            self._verify_hash = hself+1
        else:
            self._verify_hash = hself

    def reset(self,recursive=False,changedstate=False):
        self.changed = changedstate
        if recursive:
            for child in self._children:
                child.reset(recursive=recursive,changedstate=changedstate)


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
from abc import ABCMeta, abstractmethod
from future.utils import with_metaclass

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


class BaseNode(collections.MutableMapping):
    # Original idea:
    # http://code.activestate.com/recipes/577434-nested-contexts-a-chain-of-mapping-objects/

    def __init__(self):
        self._hash_map_cached = None
        self._hash_all_last = None
        self._newkeys = set()
        self._subnodectr = 0

    @property
    @abstractmethod
    def name(self):
        pass

    @property
    @abstractmethod
    def parent(self):
        pass
        
    @property
    @abstractmethod
    def children(self):
        pass
        
    @property
    @abstractmethod
    def _prefix(self):
        pass

    @property
    @abstractmethod
    def _local_map(self):
        pass
        
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
        """Update map locally and remove children override of the keys
        """
        d = self.dtype(*args, **kwargs)
        self.update(d)
        for key in d:
            self.delete_recursive(key,local=False)

    def _annotate_key(self,key,localcontext=True):
        return AnnotatedKey(key,localcontext,changed=key in self._newkeys)
    
    def _map(self,annotated=False,localcontext=True):
        _map = self._local_map
        if annotated:
            if self.dtype is dict:
                return {self._annotate_key(k,localcontext=localcontext):v for k,v in _map.items()}
            else:
                return self.dtype(map(lambda item:(self._annotate_key(item[0],localcontext=localcontext),item[1]),_map.items()))
        else:
            return _map 

    def _maps(self,withself=True,annotated=False):
        """Generator of all node maps in ascending order of inheritance (self first)
        """
        if withself:
            yield self._map(annotated=annotated)
        node = self.parent
        while node is not None:
            yield node._map(annotated=annotated,localcontext=False)
            node = node.parent

    def iter_up(self,withself=True):
        """Generator of all nodes in ascending order (self first)
        """
        if withself:
            yield self
        node = self.parent
        while node is not None:
            yield node
            node = node.parent

    def iter_down(self,withself=True):
        if withself:
            yield self
        for child in self.children:
            for node in child.iter_down():
                yield node
    
    @property
    def root(self):
        for node in self.iter_up():
            pass
        return node

    def _getitem(self, key, withself=True):
        m = {}
        for m in self._maps(withself=withself):
            if key in m:
                break
        return m[key]
        
    def _getitem_local(self, key):
        return self._map()[key]
    
    def __getitem__(self, key):
        return self._getitem(key)

    @contextlib.contextmanager
    def _changekey(self,key):
        try:
            yield self._map()
        except Exception as e:
            raise e
        else:
            self._onupdate_key(key)

    def _onupdate_key(self,key):
        self._newkeys.add(key)
        self.rehash()

    def _ensure_local_subnodes(self):
        cmap = self._map()
        for k,v in self.items(withself=False):
            if isinstance(v,SubNode) and k not in cmap:
                self._local_subnode(cmap,k)

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

        if isinstance(value,collections.MutableMapping) and not isinstance(value,Node):
            self._setitem_subnode(key,value)
        else:
            self._setitem(key,value)
    
    def _setitem(self,key,value):
        with self._changekey(key) as cmap:
            cmap[key] = value
                
    def _setitem_subnode(self,key,value):
        if isinstance(value,BaseSubNode):
            value = value.todict(full=False)
        subnode = self._getitem_subnode(key)
        subnode.update(value)
        if subnode.changed:
            self._onupdate_key(key)
    
    def _getitem_subnode(self,key):
        try:
            subnode = self._getitem_local(key)
        except KeyError:
            subnode = SubNode(key,self)
            with self._changekey(key) as cmap:
                cmap[key] = subnode
        else:
            if not isinstance(subnode,BaseSubNode):
                raise RuntimeError("Key {} of {} should be undefined or a SubNode".format(key,self.name))
        return subnode
        
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
        return self.iter()

    def iter(self,withself=True,annotated=False):
        # Iterate over keys
        it = itertools.chain.from_iterable(self._maps(withself=withself,annotated=annotated))
        if annotated:
            key = lambda key:key.original
        else:
            key = None
        return listtools.unique_everseen(it,key=key)
        
    def keys(self,**kwargs):
        return self.iter(**kwargs)
    
    def items(self,withself=True,annotated=False):
        if annotated:
            return map(lambda key:(key,self._getitem(key.original,withself=withself)),self.keys(withself=withself,annotated=annotated))
        else:
            if withself:
                return super(BaseNode,self).items()
            else:
                return map(lambda key:(key,self._getitem(key,withself=withself)),self.keys(withself=withself,annotated=annotated))
                
    def __contains__(self, key):
        return any(key in m for m in self._maps())

    def isinheritedkey(self, key):
        return any(key in m for m in self._maps()[1:])

    def islocalkey(self, key):
        return key in self._map()

    def tree(self,level=0,onlychanged=False):
        """Show inheritance tree of this node
        """
        tab = "  "*level
        if tab:
            tab += "> "

        lst = [child.tree(level=level+1,onlychanged=onlychanged) for child in self.children]
        lst.insert(0,tab+self._prefix+self._repr_local(onlychanged=onlychanged))
        
        return "\n".join(lst)
    
    def _item_changed(self,k,v):
        if k in self._newkeys:
            return True
        if isinstance(v,BaseNode):
            return v.changed
        return False
    
    def _repr_local(self, onlychanged=False):
        if onlychanged:
            lst = ["{}: {}".format(k,v) for k,v in self._map().items() if self._item_changed(k,v)]
        else:
            lst = ["{}{}: {}".format(k,"(*)" if self._item_changed(k,v) else "",v) for k,v in self._map().items()]
            
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
        return  " -> ".join([node._repr_local() for node in self.iter_up()])

    def __repr__(self):
        return self._prefix + self._repr_all()

    @property
    def _hash_local(self):
        # disabled caching
        return hashing.calchash(self._map())
        #if not self._hash_map_cached:
        #    self._hash_map_cached = hashing.calchash(self._map())
        #return self._hash_map_cached

    def __hash__(self):
        # Hash of all local map up the tree
        return hash((self._hash_local,hash(self.parent)))
  
    @property
    def changed(self):
        return hash(self)!=self._hash_all_last

    @changed.setter
    def changed(self,changed):
        if changed:
            self._hash_all_last = None
            self._newkeys = set(self._map())
        else:
            self._hash_all_last = hash(self)
            self._newkeys = set()
        
    def reset_state(self,changed=False,recursive=False):
        self.changed = changed
        for v in self.values():
            if isinstance(v,SubNode):
                v.reset_state(changed=changed)
        if recursive:
            for child in self.children:
                child.reset_state(changed=changed,recursive=recursive)

class Node(BaseNode):

    def __init__(self, parent=None, name=None):
        super(Node,self).__init__()

        # Local parameters
        self._name = name
        self._data = {}

        # Link to other nodes
        self._parent = None
        self._children = []
        self.parent = parent

    @property
    def name(self):
        if self._name:
            return self._name
        else:
            return str(id(self))
    
    @property
    def _prefix(self):
        return "{} = ".format(self.name)
    
    @property
    def _local_map(self):
        return self._data
    
    @property
    def children(self):
        return self._children
    
    @children.setter
    def children(self,lst):
        self._children = lst
    
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
            if id(node) in map(id,self.iter_down()):
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
        self.parent = None
        return newparent

class BaseSubNode(BaseNode):

    def __init__(self):
        super(BaseSubNode,self).__init__()

    @property
    def _prefix(self):
        return ""
        
class SubNode(BaseSubNode):

    def __init__(self, key, parent):
        self._key = key
        self._data = {}
        self._direct_parent = parent
        super(SubNode,self).__init__()

    @property
    def name(self):
        return self._key
   
    @property
    def _reference(self):
        return SubNodeReferenceUp(self)
        
    @property
    def parent(self):
        return self._reference.parent
        
    @property
    def children(self):
        return self._reference.children
    
    @property
    def _local_map(self):
        return self._data
    
class SubNodeReference(BaseSubNode):

    def __init__(self):
        super(SubNodeReference,self).__init__()
        
    @property
    def parent(self):
        if self.node is None:
            return None
        else:
            return SubNodeReferenceDown(self.node.parent,self.subkeys)
    
    @property
    def children(self):
        if self.node is not None:
            for child in self.node.children:
                yield SubNodeReferenceDown(child,self.subkeys)
                
    @property
    def _local_map(self):
        subnode = self.subnode
        if subnode is None:
            return {}
        else:
            return subnode._local_map
    
    @property
    @abstractmethod
    def subnode(self):
        pass
    
    @property
    @abstractmethod
    def node(self):
        pass
        
    @property
    @abstractmethod
    def subkeys(self):
        pass

class SubNodeReferenceUp(SubNodeReference):

    def __init__(self,subnode):
        self._subnode = subnode
        self._getnode()
        super(SubNodeReferenceUp,self).__init__()

    def _getnode(self):
        node = self._subnode
        subkeys = []
        while isinstance(node,SubNode):
            subkeys.append(node._key)
            node = node._direct_parent
        self._node = node
        self._subkeys = reversed(subkeys)
        
    @property
    def subnode(self):
        return self._subnode
    
    @property
    def node(self):
        return self._node
        
    @property
    def subkeys(self):
        return self._subkeys
        
class SubNodeReferenceDown(SubNodeReference):

    def __init__(self,node,subkeys):
        self._node = node
        self._subkeys = subkeys
        super(SubNodeReferenceDown,self).__init__()

    @property
    def subnode(self):
        if self.node is None:
            return None
        subnode = self.node
        for key in self.subkeys:
            subnode = subnode._getitem_subnode(key)
        return subnode

    @property
    def node(self):
        return self._node
        
    @property
    def subkeys(self):
        return self._subkeys
        

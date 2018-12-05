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

from copy import deepcopy
from abc import ABCMeta,abstractmethod
from future.utils import with_metaclass
import logging
from contextlib import contextmanager

from ..utils import instance
from ..utils import timing
from ..utils import hashing
from ..io.utils import randomstring

logger = logging.getLogger(__name__)

class TaskException(Exception):
    pass

class MissingParameter(TaskException):
    pass

class ParameterError(TaskException):
    pass

class Task(with_metaclass(ABCMeta,object)):
    """Task who's output is a single Nexus HDF5 path
    """
    
    def __init__(self,dependencies=None,nxparent=None,**parameters):
        """
        Args:
            dependencies(Optional(Task|Path))
            nxparent(Optional(Path))
        """
        self._tempname = randomstring()
        self.parameters = parameters
        self.dependencies = dependencies
        self._nxparent = nxparent
        if not self.dependencies and not nxparent:
            raise ValueError('Specify "nxparent" when task does not have dependencies')
    
    def __str__(self):
        return "Task '{}'".format(self.name)
    
    def run(self):
        """Creates the output atomically
        """
        with timing.timeit_logger(logger,name=str(self)):
            if self.done:
                logger.info('{} already done'.format(self))
            else:
                logger.info('{} started ...'.format(self))
                
                # Make sure dependencies are done
                for task in self.dependencies:
                    if not task.done:
                        return
                
                # Process and create result atomically
                with self._atomic_context():
                    self._execute()

    @property
    def output(self):
        return self.nxparent[self.name]
    
    @property
    def done(self):
        """A task is done when
            - the dependencies exist
            - the output exists
        """
        for output in self.previous_outputs:
            if not output.exists:
                return False
        return self.output.exists
    
    @property
    def checksum(self):
        hashes = [dependency.checksum for dependency in self.dependencies]
        hashes.append(hashing.calcjhash(self.parameters))
        return hashing.mergejhash(*hashes)
        
    @property
    def parameters(self):
        return self._parameters
    
    @parameters.setter
    def parameters(self,value):
        self._parameters = deepcopy(value)
        self._parameters_defaults()
        lst = self._parameters_filter()
        if lst:
            self._parameters = {k:v for k,v in self._parameters.items() if k in lst}
    
    @property
    def dependencies(self):
        return self._dependencies
        
    @dependencies.setter
    def dependencies(self,value):
        if instance.isarray(value):
            self._dependencies = value
        else:
            if value is None:
                self._dependencies = []
            else:
                self._dependencies = [value]
                
    @property
    def previous_outputs(self):
        ret = []
        for dependency in self.dependencies:
            if isinstance(dependency,Task):
                ret.append(dependency.output)
            else:
                ret.append(dependency)
        return ret

    def _parameters_defaults(self):
        self._required_parameters('method')
        self.parameters['name'] = self.parameters.get('name',self.method)

    def _required_parameters(self,*params):
        parameters = self.parameters
        for p in params:
            if p not in parameters:
                raise MissingParameter(p)
    
    def _parameters_filter(self):
        return ['name','method']
    
    @property
    def name(self):
        return self.parameters['name']
        
    @property
    def method(self):
        return self.parameters.get('method',None)

    @property
    def nxparent(self):
        if self.dependencies:
            return self.previous_outputs[-1].parent
        else:
            return self._nxparent
        
    @abstractmethod
    def _execute(self):
        pass
    
    @abstractmethod
    @contextmanager
    def _atomic_context(self):
        pass
        
    
def task(**parameters):
    method = parameters.get('method',None)
    if method=='crop':
        from .nxcrop import Task
    elif method=='replace':
        from .nxreplace import Task
    elif method=='minlog':
        from .nxminlog import Task
    elif method=='align':
        from .nxalign import Task
    elif method=='expression':
        from .nxexpression import Task
    elif method=='resample':
        from .nxresample import Task
    elif method=='xrf':
        from .nxxrf import Task
    elif method=='fullfield':
        from .nxfullfield import Task
    elif method=='xiaedftonx':
        from .nxxiaedf import Task
    else:
        path = parameters.get('path',None)
        if path is None:
            raise TaskException('Unknown task {}'.format(repr(method)))
        if path.is_nxclass('NXprocess'):
            from .nxprocesswrap import Task
        else:
            from .nxwrap import Task
    return Task(**parameters)


def nxpathtotask(path):
    if path.is_nxclass('NXprocess'):
        parameters = path.config.read()
        if 'method' not in parameters and 'name' not in parameters:
            parameters['path'] = path
    else:
        parameters = {'path':path}
    nxparent = path.parent
    dependencies = [path for path in path.dependencies]
    return task(dependencies=dependencies,nxparent=nxparent,**parameters)
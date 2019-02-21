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
import traceback

from ..utils import instance
from ..utils import timing
from ..utils import hashing
from ..utils.signalhandling import DelaySignalsContext

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
    
    def __init__(self,dependencies=None,outputparent=None,**parameters):
        """
        Args:
            dependencies(Optional(Task|Path))
            outputparent(Optional(Path))
        """
        self.outputparent = outputparent
        self.dependencies = dependencies
        self.parameters = parameters
        if not self.hasdependencies and self.outputparent is None:
            raise ValueError('Specify "outputparent" when task does not have dependencies')
        
    def __str__(self):
        return "Task '{}'".format(self.output.name)
    
    def __repr__(self):
        return self.__str__()
    
    def run(self):
        """Creates the output atomically
        """
        if self.dependencies_done:
            self._ensure_outputdeviceparent()
            if self.exists:
                logger.info('{} already done'.format(self))
            else:
                logger.info('{} started ...'.format(self))
                with timing.timeit_logger(logger,name=str(self)):
                    with self._atomic_context():
                        self._execute()
        else:
            logger.warning('{} not executed (missing dependencies)'.format(self))

    def _atomic_context(self):
        return DelaySignalsContext(setup=self._atomic_context_enter,
                                   teardown=self._atomic_context_exit)

    def _atomic_context_enter(self):
        pass

    def _atomic_context_exit(self, exc_type, exc_value, exc_traceback):
        if exc_type:
            logger.error(''.join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
        return 1 # Exception is handled (do not raise it)

    @property
    def output(self):
        return self.outputparent[self.outputname]
    
    @property
    def dependencies_done(self):
        for dependency in self.dependencies:
            if hasattr(dependency,'output'):
                if not dependency.done:
                    return False
            else:
                if not dependency.exists:
                    return False
        return True
    
    @property
    def exists(self):
        return self.output.exists
            
    @property
    def done(self):
        return self.dependencies_done and self.exists
        
    @property
    def checksum(self):
        hashes = [dependency.checksum for dependency in self.dependencies]
        hashes.append(hashing.calchash(self.parameters))
        return hashing.mergehash(*hashes)
        
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
        self.default_outputparent = None
    
    @property
    def default_outputparent(self):
        if self._default_outputparent is None:
            try:
                previous = self.previous_outputs
            except AttributeError:
                pass
            else:
                if previous:
                    self._default_outputparent = previous[-1].parent
                else:
                    self._default_outputparent = None
        return self._default_outputparent
        
    @default_outputparent.setter
    def default_outputparent(self,value):
        self._default_outputparent = value
    
    @property
    def previous_outputs(self):
        ret = []
        for dependency in self.dependencies:
            if isinstance(dependency,Task):
                ret.append(dependency.output)
            else:
                ret.append(dependency)
        return ret

    @property
    def hasdependencies(self):
        return bool([x for x in self.dependencies])

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
    def method(self):
        return self.parameters['method']

    @property
    def outputname(self):
        return self.parameters['name']

    @property
    def outputparent(self):
        if self._outputparent is None:
            return self.default_outputparent
        else:
            return self._outputparent
    
    @outputparent.setter
    def outputparent(self,value):
        self._outputparent = value

    def _ensure_outputdeviceparent(self):
        device = self.outputparent.device
        if device:
            device.parent.mkdir()
        else:
            self.outputparent.mkdir()
    
    @abstractmethod
    def _execute(self):
        pass

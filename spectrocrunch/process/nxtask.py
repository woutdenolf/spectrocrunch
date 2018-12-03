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
import traceback
import logging
import contextlib

from . import nxutils
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
    
    def __init__(self,dependencies=None,nxentry=None,**parameters):
        """
        Args:
            dependencies(Optional(Task|Path))
        """
        self._tempname = randomstring()
        self.parameters = parameters
        self.dependencies = dependencies
        self._nxentry = nxentry
        self.nxprocess = None
        self.nxresults = None
        if not self.dependencies and not nxentry:
            raise ValueError('Specify "nxentry" when task does not have dependencies')
            
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
    def checksum(self):
        hashes = [dependency.checksum for dependency in self.dependencies]
        hashes.append(hashing.calcjhash(self.parameters))
        return hashing.mergejhash(*hashes)
    
    @property
    def dependencies(self):
        return self._dependencies
    
    @property
    def previous_outputs(self):
        ret = []
        for dependency in self.dependencies:
            if isinstance(dependency,Task):
                ret.append(dependency.output)
            else:
                ret.append(dependency)
        return ret
        
    @dependencies.setter
    def dependencies(self,value):
        if instance.isarray(value):
            self._dependencies = value
        else:
            if value is None:
                self._dependencies = []
            else:
                self._dependencies = [value]

    def _parameters_defaults(self):
        self._required_parameters('method')
        self.parameters['name'] = self.parameters.get('name',self.method)

    def _required_parameters(self,*params):
        parameters = self.parameters
        for p in params:
            if p not in parameters:
                raise MissingParameter(p)
    
    def _parameters_filter(self):
        return ['name','method','default']
    
    @property
    def name(self):
        return self.parameters['name']
        
    @property
    def method(self):
        return self.parameters.get('method',None)
    
    @property
    def default(self):
        return self.parameters.get('default',None)
    
    def __str__(self):
        return "Task '{}'".format(self.name)
    
    def run(self):
        """Creates an NXprocess group atomically
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
                with self._atomic_nxprocess():
                    self._execute()

    @contextlib.contextmanager
    def _atomic_nxprocess(self):
        """This is atomic if h5py.Group.move is atomic
        """
        self.nxprocess,_ = self.nxentry.nxprocess(self._tempname,
                                    parameters=self.parameters,
                                    dependencies=list(self.previous_outputs))
        self.nxresults = self.nxprocess.results
        try:
            yield
        except Exception:
            logger.error(traceback.format_exc())
            self.nxprocess.remove(recursive=True)
        else:
            nxprocess = self.nxprocess.rename(self.output)
            if self.default:
                nxutils.set_default(nxprocess,self.default)
            nxprocess.updated()
        finally:
            self.nxprocess = None
            self.nxresults = None
            
    @property
    def nxentry(self):
        if self.dependencies:
            return self.previous_outputs[-1].nxentry()
        else:
            return self._nxentry

    @property
    def output(self):
        return self.nxentry[self.name]
    
    @property
    def done(self):
        """A task is done when
            - the dependencies exist
            - the output exists with the same name and parameters
        """
        dependencies = []
        for output in self.previous_outputs:
            if output.exists:
                dependencies.append(output)
            else:
                return False
        _,exists = self.nxentry.nxprocess_exists(self.name,
                                    parameters=self.parameters,
                                    dependencies=dependencies)
        return exists
    
    @abstractmethod
    def _execute(self):
        pass

def newtask(**parameters):
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
    else:
        from .nxwrap import Task
    return Task(**parameters)

def nxprocesstotask(nxprocess):
    parameters = nxprocess.config.read()
    if 'method' not in parameters and 'name' not in parameters:
        parameters['nxprocess'] = nxprocess
    nxentry = nxprocess.nxentry()
    dependencies = [path for path in nxprocess.dependencies]
    return newtask(dependencies=dependencies,nxentry=nxentry,**parameters)

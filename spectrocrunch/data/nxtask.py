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
from ..io.utils import randomstring

logger = logging.getLogger(__name__)

class TaskException(Exception):
    pass

class MissingParameter(TaskException):
    pass

class ParameterError(TaskException):
    pass
    
class Task(with_metaclass(ABCMeta,object)):
    
    def __init__(self,previous=None,nxentry=None,**parameters):
        self._tempname = randomstring()
        self.parameters = parameters
        self.previous = previous
        self._nxentry = nxentry
        self.nxprocess = None
        self.nxresults = None
        if not self.previous and not nxentry:
            raise ValueError('Specify "nxentry" when task is not based on a previous NXprocess')
            
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
    def previous(self):
        return self._previous
    
    @previous.setter
    def previous(self,value):
        if instance.isarray(value):
            self._previous = value
        else:
            if value is None:
                self._previous = []
            else:
                self._previous = [value]
                
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
            fmt = "Task {} {{}}".format(self)
            if self.done:
                logger.info(fmt.format("already done"))
            else:
                logger.info(fmt.format("started ..."))
                
                # Make sure previous tasks are executed
                for nxprocess in self.previous:
                    if not nxprocess.exists:
                        return
                
                with self._atomic_nxprocess():
                    self._execute()

    @contextlib.contextmanager
    def _atomic_nxprocess(self):
        """This is atomic if h5py.Group.move is atomic
        """
        self.nxprocess,_ = self.nxentry.nxprocess(self._tempname,
                                    parameters=self.parameters,
                                    previous=self.previous)
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
        if self.previous:
            return self.previous[-1].nxentry()
        else:
            return self._nxentry

    @property
    def output(self):
        return self.nxentry[self.name]
    
    @property
    def done(self):
        """A task is done when the output exists with the same name and parameters
        """
        _,exists = self.nxentry.nxprocess_exists(self.name,
                                    parameters=self.parameters,
                                    previous=self.previous)
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
    else:
        raise ParameterError('Unknown method {}'.format(repr(method)))
    return Task(**parameters)

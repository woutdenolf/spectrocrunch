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

from . import nxutils

class TaskException(Exception):
    pass

class MissingParameter(TaskException):
    pass

class ParameterError(TaskException):
    pass
    
class Task(with_metaclass(ABCMeta,object)):
    
    def __init__(self,parameters,previous):
        self.parameters = parameters
        self.previous = previous
        
    @property
    def parameters(self):
        return self._parameters
    
    @parameters.setter
    def parameters(self,value):
        self._parameters = deepcopy(value)
        self._parameters_defaults()
        lst = self._parameters_filter()
        self._parameters = {k:v for k,v in self._parameters.items() if k in lst}
    
    def _parameters_defaults(self):
        self.parameters["name"] = self.parameters.get('name',self.method)

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
    
    def run(self):
        nxprocess = self._execute()
        if nxprocess and self.default:
            nxutils.set_default(nxprocess,self.default)
        return nxprocess
    
    @abstractmethod
    def _execute(self):
        pass

def newtask(parameters,previous):
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
    else:
        raise ParameterError('Unknown method {}'.format(repr(method)))
    return Task(parameters,previous)

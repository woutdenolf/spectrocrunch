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

from . import nxtask

class Task(nxtask.Task):
    """Used to wrap an NXprocess which is not produced by an NXtask
    """
    
    def __init__(self,nxprocess=None,**kwargs):
        super(Task,self).__init__(**kwargs)
        if 'name' in self.parameters:
            self.nxprocess = self.nxentry[self.parameters['name']]
            self.nxresults = self.nxprocess.results
        else:
            if nxprocess is None:
                raise ValueError('Provide "nxprocess" to the wrapper task')
            self.nxprocess = nxprocess
            self.nxresults = nxprocess.results

    def _parameters_defaults(self):
        pass
    
    def _parameters_filter(self):
        pass
        
    def _execute(self):
        pass

    @property
    def name(self):
        return self.nxprocess.name

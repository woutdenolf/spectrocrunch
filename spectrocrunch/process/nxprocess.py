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

import logging
from contextlib import contextmanager
import traceback

from . import nxutils
from . import nxtask

logger = logging.getLogger(__name__)

class Task(nxtask.Task):
    """Task who's output is a single NXprocess
    """
    
    def __init__(self,**kwargs):
        super(Task,self).__init__(**kwargs)
        self.nxprocess = None

    @property
    def output(self):
        return self.nxentry[self.name]
    
    @property
    def done(self):
        """A task is done when
            - the dependencies exist
            - the output exists
            - the parameters have the expected hash
        """
        if not super(Task,self).done:
            return False
        _,exists = self.nxentry.nxprocess_exists(self.name,
                                    parameters=self.parameters,
                                    dependencies=self.dependencies)
        return exists
        
    def _parameters_filter(self):
        return super(Task,self)._parameters_filter()+['default']

    @property
    def default(self):
        return self.parameters.get('default',None)

    @contextmanager
    def _atomic_context(self):
        """This is atomic if h5py.Group.move is atomic
        """
        self.nxprocess,_ = self.nxentry.nxprocess(self._tempname,
                                    parameters=self.parameters,
                                    dependencies=list(self.previous_outputs))
        try:
            yield
        except Exception:
            logger.error(traceback.format_exc())
            self.nxprocess.remove(recursive=True)
        else:
            self.nxprocess = self.nxprocess.rename(self.output)
            if self.default:
                nxutils.set_default(self.nxprocess,self.default)
            self.nxprocess.updated()
        finally:
            self.nxprocess = None
    
    @property
    def nxentry(self):
        return self.nxparent.nxentry()
            
    @property
    def nxresults(self):
        if self.nxprocess is None:
            return None
        else:
            return self.nxprocess.results

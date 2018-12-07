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
import traceback

from . import nxutils
from . import basetask
from ..io import nxfs

logger = logging.getLogger(__name__)


class Task(basetask.Task):
    """Task who's output is a single NXprocess
    """
    
    def __init__(self,**kwargs):
        super(Task,self).__init__(**kwargs)
        self.nxprocess = None

    @property
    def exists(self):
        _,exists = self.outputparent.nxprocess_exists(self.name,
                                        parameters=self.parameters,
                                        dependencies=self.dependencies)
        return exists
    
    @property
    def _run_alreadydone(self):
        try:
            return self.exists
        except nxfs.NexusProcessWrongHash:
            self.outputparent = self.outputparent.new_nxentry()
            return self.exists
        
    def _parameters_filter(self):
        return super(Task,self)._parameters_filter()+['default']

    @property
    def default(self):
        return self.parameters.get('default',None)

    def _atomic_context_enter(self):
        """This is atomic if h5py.Group.move is atomic
        """
        self.nxprocess,_ = self.outputparent.nxprocess(self._tempname,
                                    parameters=self.parameters,
                                    dependencies=list(self.previous_outputs))

    def _atomic_context_exit(self, exc_type, exc_value, exc_traceback):
        if exc_type:
            logger.error(''.join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
            self.nxprocess.remove(recursive=True)
        else:
            self.nxprocess = self.nxprocess.rename(self.output)
            if self.default:
                nxutils.set_default(self.nxprocess,self.default)
            self.nxprocess.updated()
        self.nxprocess = None
        return 1 # Exception is handled (do not raise it)

    @property
    def nxresults(self):
        if self.nxprocess is None:
            return None
        else:
            return self.nxprocess.results

    def _ensure_outputparent(self):
        super(Task,self)._ensure_outputparent()
        outputparent = self.outputparent
        self.outputparent = outputparent.nxentry(name=outputparent.name)

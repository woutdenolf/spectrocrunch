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
from . import target
from ..io import nxfs
from ..io import fs
from ..io.utils import randomstring
from ..utils import incremental_naming

logger = logging.getLogger(__name__)

class Task(basetask.Task):
    """Task who's output is a single self.temp_nxprocess
    """

    def __init__(self,**kwargs):
        super(Task,self).__init__(**kwargs)
        self._init_outputname()

    def _init_outputname(self):
        self._outputname = incremental_naming.Name(self.parameters['name'])
        
    def _update_output_name(self):
        entry = self.outputparent
        process = entry.find_nxprocess(name=self.outputname,parameters=self.parameters,
                                       dependencies=self.dependencies,searchallentries=False)
        if process is None:
            while entry[self.outputname].exists:
                self._outputname += 1
        else:
            self._outputname = incremental_naming.Name(process.name)
            
    @property
    def output(self):
        self._update_output_name()
        return self.outputparent[self.outputname]
    
    @property
    def outputname(self):
        return str(self._outputname)
        
    def _parameters_filter(self):
        return super(Task,self)._parameters_filter()+['default']

    @property
    def default(self):
        return self.parameters.get('default',None)

    def _atomic_context_enter(self):
        while True:
            try:
                self.temp_nxprocess = self.outputparent.nxprocess(randomstring(),noincrement=True,
                                                                parameters=self.parameters,
                                                                dependencies=list(self.previous_outputs))
            except fs.AlreadyExists:
                pass # already exists
            else:
                break
        
    @property
    def temp_nxresults(self):
        return self.temp_nxprocess.results
        
    @property
    def temp_name(self):
        return self.temp_nxprocess.name
        
    def _atomic_context_exit(self, exc_type, exc_value, exc_traceback):
        if exc_type:
            logger.error(''.join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
            self.removeoutput()
        else:
            self.renameoutput()
            nxutils.set_default(self.output,self.default)
        self.temp_nxprocess = None
        return 1 # Exception is handled (do not raise it)

    def removeoutput(self):
        self.temp_nxprocess.remove(recursive=True)
    
    def renameoutput(self):
        self._init_outputname()
        while self.temp_nxprocess.exists:
            try:
                self.temp_nxprocess.rename(self.output)
            except fs.AlreadyExists:
                if self.output.exists:
                    # Already done by someone else
                    self.removeoutput()

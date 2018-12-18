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
from ..io import fs

logger = logging.getLogger(__name__)

class Task(basetask.Task):
    """Task who's output is a single self.temp_nxprocess
    """

    def __init__(self,**kwargs):
        super(Task,self).__init__(**kwargs)
        self.temp_nxprocess = None
        self._outputcounter = 0
    
    @property
    def output(self):
        return self.outputparent.get_nxprocess(self.outputname,parameters=self.parameters,
                                                dependencies=self.dependencies,allentries=False)
    
    @property
    def outputname(self):
        # Add a counter when needed:
        #   'align' -> 'align.1'
        #   'align0001' -> 'align0001'
        name = self.parameters['name']
        fmt,num,nonumber = nxfs.nxprocess_fmtname(name)
        if nonumber:
            num += 1
        num +=self._outputcounter
        return fmt.format(num)

    def _parameters_filter(self):
        return super(Task,self)._parameters_filter()+['default']

    @property
    def default(self):
        return self.parameters.get('default',None)

    def _atomic_context_enter(self):
        """This is atomic if h5py.Group.move is atomic
        """
        self._outputcounter = 0
        self.temp_nxprocess = self.output.parent.nxprocess(self._tempname,
                                                        parameters=self.parameters,
                                                        dependencies=list(self.previous_outputs))
        
    def _atomic_context_exit(self, exc_type, exc_value, exc_traceback):
        if exc_type:
            logger.error(''.join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
            self.temp_nxprocess.remove(recursive=True)
        else:
            while self.temp_nxprocess.exists:
                checksum = self.temp_nxprocess.checksum
                try:
                    self.temp_nxprocess.rename(self.output)
                except fs.AlreadyExists:
                    if self.output.checksum==checksum:
                        self.temp_nxprocess.remove(recursive=True)
                    else:
                        self._outputcounter += 1
            if self.default:
                nxutils.set_default(self.output,self.default)
            else:
                self.output.mark_default()
            self.output.updated()
        self.temp_nxprocess = None
        return 1 # Exception is handled (do not raise it)

    @property
    def temp_nxresults(self):
        if self.temp_nxprocess is None:
            return None
        else:
            return self.temp_nxprocess.results

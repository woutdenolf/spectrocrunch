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
from ..io import target
from ..io import nxfs
from ..io import fs
from ..io.utils import randomstring

logger = logging.getLogger(__name__)


class Task(basetask.Task):
    """Task who's output is a single NXprocess
    """

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.optional_parameters.add('default')

    @property
    def default(self):
        return self.parameters.get('default', None)

    def _atomic_context_enter(self):
        while True:
            try:
                self.temp_nxprocess = self.outputparent.nxprocess(randomstring(),
                                                                  noincrement=True,
                                                                  parameters=self.parameters,
                                                                  dependencies=list(self.previous_outputs))
            except fs.AlreadyExists:
                pass  # already exists
            else:
                break

    @property
    def outputname(self):
        """
        NXprocess name (non-existent unless locked)
        """
        if not hasattr(self, '_target'):
            self._target = target.Name(self.parameters['name'])
        if self._target.locked:
            return str(self._target)
        else:
            self._target.reset()
        entry = self.outputparent
        if not entry.exists:
            return str(self._target)
        process = entry.find_nxprocess(name=str(self._target),
                                       parameters=self.parameters,
                                       dependencies=self.dependencies,
                                       searchallentries=False)
        if process is None:
            while entry[str(self._target)].exists:
                self._target += 1
        else:
            self._target = target.Name(process.name)
        return str(self._target)

    def _lock_outputname(self):
        self._target.lock()

    @property
    def temp_nxresults(self):
        return self.temp_nxprocess.results

    @property
    def temp_outputname(self):
        return self.temp_nxprocess.name

    @property
    def localpath(self):
        return self.outputparent.device.parent

    @property
    def temp_localpath(self):
        return self.localpath[self.temp_outputname]

    @property
    def output_localpath(self):
        name = self.outputparent.name + '_' + self.outputname
        return self.localpath[name]

    def _atomic_context_exit(self, exc_type, exc_value, exc_traceback):
        if exc_type:
            logger.error(''.join(traceback.format_exception(
                exc_type, exc_value, exc_traceback)))
            self._remove_temp_output()
        else:
            self._rename_temp_output()
            nxutils.set_default(self.output, self.default)
        self.temp_nxprocess = None
        return 1  # Exception is handled (do not raise it)

    def _remove_temp_output(self):
        self.temp_nxprocess.remove(recursive=True)
        self.temp_localpath.remove(recursive=True)

    def _rename_temp_output(self):
        while self.temp_nxprocess.exists:
            try:
                old, new = self.temp_nxprocess, self.output
                old.rename(new)
            except fs.AlreadyExists:
                if self.output.exists:
                    # Already done by someone else
                    self._remove_temp_output()
            else:
                logger.info('Renamed {} to {}'.format(old, new))
                self._lock_outputname()
        old = self.temp_localpath
        if old.exists:
            new = self.output_localpath
            old.move(new, force=True)
            logger.info('Renamed {} to {}'.format(old, new))

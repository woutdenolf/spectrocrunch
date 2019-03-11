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

import re
import logging
import traceback

from . import nxqxrf_dependent
from ..io import xiaedf
from ..io import xiaedftonexus
from ..io import fs
from ..io.utils import randomstring

logger = logging.getLogger(__name__)


class Task(nxqxrf_dependent.Task):
    """Converts XIA edf to an NXentry
    """

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.required_parameters |= {
            'path',
            'radix',
            'number',
            'instrument'
        }

        self.optional_parameters |= {
            'include_counters',
            'exclude_counters',
            'fluxid',
            'transmissionid'
        }

    def _atomic_context_enter(self):
        name = randomstring()
        root = self.outputparent
        while root[name].exists:
            name = randomstring()
        self.temp_nxentry = root.nxentry(name=name)

    def _atomic_context_exit(self, exc_type, exc_value, exc_traceback):
        if exc_type:
            logger.error(''.join(traceback.format_exception(
                exc_type, exc_value, exc_traceback)))
            self.temp_nxentry.remove(recursive=True)
        else:
            self.temp_nxentry = self.temp_nxentry.renameremove(self.output)
            self.temp_nxentry.mark_default()
        self.temp_nxentry = None
        return 1  # Exception is handled (do not raise it)

    def _execute(self):
        parameters = self.parameters
        include_counters = [self._rematch_func(
            redict) for redict in parameters.get('include_counters', [])]
        exclude_counters = [self._rematch_func(
            redict) for redict in parameters.get('exclude_counters', [])]
        path, radix, number = (parameters['path'],
                               parameters['radix'],
                               parameters['number'])
        instrument = parameters.get('instrument', None)
        fluxid = parameters.get('fluxid', None)
        transmissionid = parameters.get('transmissionid', None)
        xiaimage = xiaedf.xiaimage_number(path, radix, number)
        converter = xiaedftonexus.Converter(
            nxentry=self.temp_nxentry,
            qxrfgeometry=self.qxrfgeometry,
            include_counters=include_counters,
            exclude_counters=exclude_counters,
            instrument=instrument,
            fluxid=fluxid,
            transmissionid=transmissionid)
        nxentry = converter(xiaimage)

    @property
    def name(self):
        parameters = self.parameters
        radix, number = parameters['radix'], parameters['number']
        return radix+'.{}'.format(number)

    @staticmethod
    def _rematch_func(redict):
        method = redict.get('method', 'regex')
        if method == 'regex':
            return lambda ctrname: re.match(redict['pattern'], ctrname)
        elif method == 'equal':
            return lambda ctrname: redict['value'] == ctrname
        else:
            return lambda ctrname: False

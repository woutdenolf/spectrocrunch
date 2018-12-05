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

from contextlib import contextmanager
import re
import logging
import traceback
from copy import deepcopy

from . import nxtask
from ..io import xiaedf
from ..io import xiaedftonexus

logger = logging.getLogger(__name__)

class Task(nxtask.Task):
    """Converts XIA edf output to an NXentry
    """
  
    def _parameters_defaults(self):
        super(Task,self)._parameters_defaults()
        self._required_parameters('path','radix','number','instrument')
        
    def _parameters_filter(self):
        return []
        
    @contextmanager
    def _atomic_context(self):
        """This is atomic if h5py.Group.move is atomic
        """
        self.nxentry = self._nxparent.nxentry(name=self._tempname)
        try:
            yield
        except Exception:
            logger.error(traceback.format_exc())
            self.nxentry.remove(recursive=True)
        else:
            self.nxentry = self.nxentry.rename(self.output)
            self.nxentry.updated()
        finally:
            self.nxentry = None

    def _execute(self):
        parameters = deepcopy(self.parameters)
        parameters['include_counters'] = [self._rematch_func(redict) for redict in parameters.get('include_counters',[])]
        parameters['exclude_counters'] = [self._rematch_func(redict) for redict in parameters.get('exclude_counters',[])]
        path,radix,number = parameters['path'],parameters['radix'],parameters['number']
        xiaimage = xiaedf.xiaimage_number(path,radix,number)
        converter = xiaedftonexus.Converter(nxentry=self.nxentry,**parameters)
        nxentry = converter(xiaimage)
    
    @property
    def name(self):
        parameters = self.parameters
        radix,number = parameters['radix'],parameters['number']
        return radix+'.{}'.format(number)
    
    @staticmethod
    def _rematch_func(redict):
        method = redict.get('method','regex')
        if method=='regex':
            return lambda ctrname:re.match(redict['pattern'],ctrname)
        elif method=='equal':
            return lambda ctrname:redict['value'] == ctrname
        else:
            return lambda ctrname:False
    
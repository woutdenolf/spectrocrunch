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

import contextlib
import re
import logging
import traceback

from . import nxtask
from ..io import xiaedf
from ..io import xiaedftonexus

logger = logging.getLogger(__name__)

class Task(nxtask.Task):
    """Converts XIA edf output to an NXentry
    """

    @property
    def nxroot(self):
        return self._nxparent.nxroot()
            
    def _parameters_defaults(self):
        super(Task,self)._parameters_defaults()
        self._required_parameters('path','radix','number','instrument')
        
    def _parameters_filter(self):
        return []
        
    @contextlib.contextmanager
    def _atomic_context(self):
        """This is currently not atomic!!!
        """
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
        self.nxentry = self.nxroot[self._tempname]
        parameters = self.parameters
        path,radix,number = parameters['path'],parameters['radix'],parameters['number']
        xiaimage = xiaedf.xiaimage_number(path,radix,number)
        h5file = str(self.nxroot.device)
        nxentry = self.nxentry.name
        converter = xiaedftonexus.Converter(h5file=h5file,nxentry=nxentry,**parameters)
        nxentry = converter(xiaimage)
    
    @property
    def name(self):
        parameters = self.parameters
        radix,number = parameters['radix'],parameters['number']
        return radix+'.{}'.format(number)
        

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

from . import target

class Task(basetask.Task):
    """Create scene image
    """
    
    def __init__(self,extensions,**kwargs):
        super(Task,self).__init__(**kwargs)
        self._final_output = target.TargetLocalFs(self.outputparent,self.outputname,extensions)
        
    @property
    def output(self):
        return self._final_output
    
    def _atomic_context_enter(self):
        self._temp_output = {ext:self.outputparent[self._tempname+ext] for ext in ['.png','.json']}

    def _atomic_context_exit(self, exc_type, exc_value, exc_traceback):
        if exc_type:
            logger.error(''.join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
            for path in self._temp_output.values():
                path.remove(recursive=True)
        else:
            for ext,path in self._temp_output.items():
                path.renameremove(self._final_output[ext])
        self._temp_output = None
        return 1 # Exception is handled (do not raise it)
        
    def _parameters_defaults(self):
        super(Task,self)._parameters_defaults()
        self._required_parameters('objects')
        
    def _parameters_filter(self):
        return []
        
    def _execute(self):
        parameters = self.parameters
        createscene(parameters,self._temp_output['.png'].path)
        with self._temp_output['.json'].open(mode='w') as outfile:
            json.dump(parameters, outfile)
        

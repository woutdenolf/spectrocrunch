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

from .basetask import TaskException

def create_task(**parameters):
    method = parameters.get('method', None)
    if method == 'crop':
        from .nxcrop import Task
    elif method == 'replace':
        from .nxreplace import Task
    elif method == 'minlog':
        from .nxminlog import Task
    elif method == 'align':
        from .nxalign import Task
    elif method == 'expression':
        from .nxexpression import Task
    elif method == 'resample':
        from .nxresample import Task
    elif method == 'pymca':
        from .nxpymca import Task
    elif method == 'fullfield':
        from .nxfullfield import Task
    elif method == 'xiaedftonx':
        from .nxxiaedf import Task
    elif method == 'scenevis':
        from .scenevis import Task
    elif method == 'xrfgeometry':
        from .nxqxrf import Task
    else:
        Task = parameters.pop('_task', None)
        if Task is None:
            raise TaskException(
                'Unknown task defined by parameters {}'.format(parameters))
    return Task(**parameters)


def nxpathtotask(path):
    if path.is_nxclass('NXprocess'):
        parameters = path.config.read()
        if 'method' not in parameters and 'name' not in parameters:
            parameters['name'] = path.name
        from .nxprocesswrap import Task
    else:
        from .nxwrap import Task
        parameters = {'path': path}
    outputparent = path.parent
    dependencies = [path for path in path.dependencies]
    return create_task(dependencies=dependencies, outputparent=outputparent, _task=Task, **parameters)

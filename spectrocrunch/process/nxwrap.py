# -*- coding: utf-8 -*-

from . import basetask


class Task(basetask.Task):
    """Used to wrap an HDF5 path which is not produced by nxprocess.Task"""

    def __init__(self, path=None, **kwargs):
        if path is None:
            raise basetask.TaskException('Provide "path" to the wrapper task')
        self.nxpath = path
        super(Task, self).__init__(**kwargs)

    def _parameters_defaults(self):
        pass

    def _parameters_filter(self):
        pass

    def _execute(self):
        pass

    @property
    def name(self):
        return self.nxpath.name

    @property
    def outputparent(self):
        return self.nxpath.parent

    @property
    def output(self):
        return self.nxpath

from . import nxprocess


class Task(nxprocess.Task):
    """Used to wrap an NXprocess which is produced by nxprocess.Task"""

    def _parameters_defaults(self):
        pass

    def _parameters_filter(self):
        pass

    def _execute(self):
        pass

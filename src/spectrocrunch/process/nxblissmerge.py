from . import nxentry
from .h5merge import MergedBlissMesh


class Task(nxentry.Task):
    """Merges Bliss scans (2D maps) in a 3D stack"""

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.required_parameters |= {"uris"}

    def _execute(self):
        uris = self.parameters["uris"]
        filename = self.temp_nxentry.device.path
        name = self.temp_nxentry.path
        with MergedBlissMesh(uris, filename=filename, name=name):
            pass
        self.temp_nxentry.nxprocess("merge", parameters={})

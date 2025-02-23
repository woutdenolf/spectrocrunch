import numpy as np
from . import nxregulargrid


class Task(nxregulargrid.Task):
    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.required_parameters |= {"old", "new"}

    def _process_data(self, data):
        v1 = self.parameters["old"]
        v2 = self.parameters["new"]
        if v1 is np.nan:
            data[np.isnan(data)] = v2
        else:
            data[data == v1] = v2
        return data

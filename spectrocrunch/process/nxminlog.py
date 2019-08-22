# -*- coding: utf-8 -*-

import numpy as np
from . import nxregulargrid


class Task(nxregulargrid.Task):

    def _process_data(self, data):
        return -np.log(data)

# -*- coding: utf-8 -*-
from . import nxprocess
from ..geometries import qxrf


class Task(nxprocess.Task):

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.allparams = [
            # Geometry
            'geometry',
            'init',
            # Measurements
            'fixed',
            'variable'
        ]
        self._required_parameters(*self.allparams)

    def _parameters_filter(self):
        return super(Task, self)._parameters_filter() +\
            self.allparams
        
    def _execute(self):
        parameters = self.parameters
        geom = qxrf.factory(parameters['geometry'], **parameters['init'])
        geom.batchcalibrate_diodes(parameters['fixed'], parameters['variable'])
        note = self.temp_nxresults.nxnote('geometry')
        note.write(geom)

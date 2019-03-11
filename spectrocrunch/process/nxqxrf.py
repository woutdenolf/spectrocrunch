# -*- coding: utf-8 -*-
from . import nxprocess
from ..geometries import qxrf


class Task(nxprocess.Task):

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.required_parameters |= {
            # Geometry
            'geometry',
            'init',
            # Measurements
            'fixed',
            'variable'
        }
        
    def _execute(self):
        parameters = self.parameters
        geom = qxrf.factory(parameters['geometry'], **parameters['init'])
        geom.batchcalibrate_diodes(parameters['fixed'], parameters['variable'])
        note = self.temp_nxresults.nxnote('geometry')
        note.write(geom)

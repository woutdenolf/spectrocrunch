# -*- coding: utf-8 -*-

import collections
from contextlib import ExitStack
import numpy
from . import nxmerge
from .h5merge import merge_h5groups


class Task(nxmerge.Task):
    """Tile previous NXdata's of previous tasks"""

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.required_parameters |= {"tile_shape"}

    def _execute(self):
        parameters = self.parameters
        groups = self._stack_sources()

        scan_shape = self._scan_shape(parameters["shape_parser"])
        if not scan_shape:
            scan_shape = groups["parameters"][0].shape
        tile_shape = parameters["tile_shape"]

        shape_map = {}
        axes = collections.OrderedDict()

        for name, sources in groups.items():
            for source in sources:
                shape_map[source.shape] = scan_shape
            for i, (n, mult) in enumerate(zip(scan_shape[::-1], tile_shape[::-1]), 1):
                axes["dim" + str(i)] = numpy.arange(n * mult)
            with ExitStack() as stack:
                ctx = self.temp_nxresults.open()
                dest_parent = stack.enter_context(ctx)
                nxdatas = []
                for source in sources:
                    ctx = source.open()
                    nxdata = stack.enter_context(ctx)
                    nxdatas.append(nxdata)
                merge_h5groups(
                    dest_parent,
                    name,
                    nxdatas,
                    shape_map,
                    nscandim=2,
                    tile_shape=tile_shape,
                )
            dest_parent = self.temp_nxresults[name]
            for k, v in axes.items():
                dest_parent[k].write(data=v)
            dest_parent.update_stats(axes=list(axes.keys()))

# -*- coding: utf-8 -*-

import collections
from contextlib import contextmanager, ExitStack
import numpy
from . import nxprocess
from . import basetask
from ..io.fs import Missing
from .h5merge import merge_h5groups


class Task(nxprocess.Task):
    """Stack and reshape previous NXdata's of previous tasks"""

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.optional_parameters |= {"stack_positioner", "shape"}
        parameters = self.parameters
        parameters["stack_positioner"] = parameters.get("stack_positioner", None)
        parameters["shape"] = parameters.get("shape", None)

    def _execute(self):
        parameters = self.parameters
        stack_axis_name = parameters["stack_positioner"]
        stack_axis = self._get_stack_axis(stack_axis_name)
        groups = self._stack_sources()

        idx = numpy.argsort(stack_axis)
        groups = {k: [lst[i] for i in idx] for k, lst in groups.items()}
        stack_axis = numpy.array(stack_axis)[idx]

        scan_shape = self._scan_shape()
        assert scan_shape, "Cannot extract scan shape from title"
        shape_map = {}
        axes = collections.OrderedDict()
        axes[stack_axis_name] = stack_axis

        for name, sources in groups.items():
            for source in sources:
                shape_map[source.shape] = scan_shape
            for i, n in enumerate(scan_shape[::-1], 1):
                axes["dim" + str(i)] = numpy.arange(n)
            with ExitStack() as stack:
                ctx = self.temp_nxresults.open()
                dest_parent = stack.enter_context(ctx)
                nxdatas = []
                for source in sources:
                    ctx = source.open()
                    nxdata = stack.enter_context(ctx)
                    nxdatas.append(nxdata)
                merge_h5groups(dest_parent, name, nxdatas, shape_map, nscandim=2)
            dest_parent = self.temp_nxresults[name]
            for k, v in axes.items():
                dest_parent[k].write(data=v)
            dest_parent.update_stats(axes=list(axes.keys()))

    def _get_stack_axis(self, stack_axis_name):
        pos_name = "instrument/positioners/" + stack_axis_name
        stack_axis = []
        for nxentry in self._iter_nxentry_dependencies():
            try:
                p = nxentry[pos_name].read()
            except Missing:
                p = numpy.nan
            stack_axis.append(p)
        return stack_axis

    def _stack_sources(self):
        groups = {}
        for nxprocess in self.previous_outputs:
            for nxdata in nxprocess.results.iter_is_nxclass(u"NXdata"):
                lst = groups.get(nxdata.name, None)
                if lst is None:
                    lst = groups[nxdata.name] = []
                lst.append(nxdata)
        return groups

    def _scan_shape(self):
        for nxentry in self._iter_nxentry_dependencies():
            cmd = nxentry["title"].read()
            # l2scan
            parts = cmd.split(" ")
            if parts[0] == "l2scan":
                shape = int(parts[4]), int(parts[8]) + 1
            else:
                shape = int(parts[4]) + 1, int(parts[8]) + 1
            return shape

    def _iter_nxentry_dependencies(self, dep=None):
        if dep is None:
            dep = self
        if isinstance(dep, basetask.Task):
            for dep2 in dep.dependencies:
                for dep3 in self._iter_nxentry_dependencies(dep2):
                    if not isinstance(dep3, basetask.Task):
                        if dep3.is_nxclass(u"NXentry"):
                            yield dep3
        else:
            if dep.is_nxclass(u"NXentry"):
                yield dep

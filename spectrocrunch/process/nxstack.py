# -*- coding: utf-8 -*-

import collections
from contextlib import ExitStack
import numpy
from . import nxprocess
from . import basetask
from ..io.fs import Missing
from .h5merge import merge_h5groups


class Task(nxprocess.Task):
    """Stack and reshape previous NXdata's of previous tasks"""

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.optional_parameters |= {"stack_positioner", "shape", "shape_parser"}
        parameters = self.parameters
        parameters["stack_positioner"] = parameters.get("stack_positioner", None)
        parameters["shape"] = parameters.get("shape", None)
        parameters["shape_parser"] = parameters.get("shape_parser", self._shape_parser)

    def _execute(self):
        parameters = self.parameters
        groups = self._stack_sources()

        scan_shape = self._scan_shape(parameters["shape_parser"])
        if not scan_shape:
            scan_shape = groups["parameters"][0].shape

        shape_map = {}
        axes = collections.OrderedDict()

        stack_axis_name = parameters["stack_positioner"]
        if stack_axis_name:
            stack_axis = self._get_stack_axis(stack_axis_name)
            idx = numpy.argsort(stack_axis)
            groups = {k: [lst[i] for i in idx] for k, lst in groups.items()}
            stack_axis = numpy.array(stack_axis)[idx]
            axes[stack_axis_name] = stack_axis

        for name, sources in groups.items():
            no_stack_dimension = len(sources) == 1 and not stack_axis_name
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
                merge_h5groups(
                    dest_parent,
                    name,
                    nxdatas,
                    shape_map,
                    nscandim=2,
                    no_stack_dimension=no_stack_dimension,
                )
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

    def _scan_shape(self, shape_parser):
        for nxentry in self._iter_nxentry_dependencies():
            return shape_parser(nxentry["title"].read())

    @staticmethod
    def _shape_parser(cmd):
        parts = cmd.split(" ")
        if parts[0] == "l2scan":
            shape = int(parts[4]), int(parts[8]) + 1
        else:
            try:
                shape = int(parts[4]) + 1, int(parts[8]) + 1
            except Exception:
                shape = None
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

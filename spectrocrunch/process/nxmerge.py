# -*- coding: utf-8 -*-
from . import nxprocess
from . import basetask


class Task(nxprocess.Task):
    """Merge previous NXdata's of previous tasks"""

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.optional_parameters |= {"shape", "shape_parser"}
        parameters = self.parameters
        parameters["shape"] = parameters.get("shape", None)
        parameters["shape_parser"] = parameters.get("shape_parser", self._shape_parser)

    def _stack_sources(self):
        groups = {}
        for nxprocess in self.previous_outputs:
            for nxdata in nxprocess.results.iter_is_nxclass("NXdata"):
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
        if parts[0] == "fscan2d":
            # fscan2d slowmot slowstart slowstep slownpts fastmot faststart faststep fastnpts expo1 expo2
            shape = int(parts[8]), int(parts[4])
        elif parts[0] == "l2scan":
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
                        if dep3.is_nxclass("NXentry"):
                            yield dep3
        else:
            if dep.is_nxclass("NXentry"):
                yield dep

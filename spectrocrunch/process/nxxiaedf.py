# -*- coding: utf-8 -*-

import re
import logging
import traceback

from . import nxqxrf_dependent
from ..io import xiaedf
from ..io import xiaedftonexus
from ..io import fs
from ..io.utils import randomstring

logger = logging.getLogger(__name__)


class Task(nxqxrf_dependent.Task):
    """Converts XIA edf to an NXentry"""

    def _parameters_defaults(self):
        super(Task, self)._parameters_defaults()
        self.required_parameters |= {
            "path",
            "radix",
            "number",
            "instrument",
            "include_counters",
            "exclude_counters",
            "fluxid",
            "transmissionid",
        }
        parameters = self.parameters
        parameters["fluxid"] = parameters.get("fluxid", None)
        parameters["transmissionid"] = parameters.get("transmissionid", None)
        parameters["include_counters"] = parameters.get("include_counters", [])
        parameters["exclude_counters"] = parameters.get("exclude_counters", [])

    def _atomic_context_enter(self):
        name = randomstring()
        root = self.outputparent
        while root[name].exists:
            name = randomstring()
        self.temp_nxentry = root.nxentry(name=name)

    def _atomic_context_exit(self, exc_type, exc_value, exc_traceback):
        if exc_type:
            logger.error(
                "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
            )
            self.temp_nxentry.remove(recursive=True)
        else:
            self.temp_nxentry = self.temp_nxentry.renameremove(self.output)
            self.temp_nxentry.mark_default()
        self.temp_nxentry = None
        return 1  # Exception is handled (do not raise it)

    def _execute(self):
        parameters = self.parameters
        include_counters = [
            self._rematch_func(redict) for redict in parameters["include_counters"]
        ]
        exclude_counters = [
            self._rematch_func(redict) for redict in parameters["exclude_counters"]
        ]
        path, radix, number = (
            parameters["path"],
            parameters["radix"],
            parameters["number"],
        )
        instrument = parameters["instrument"]
        fluxid = parameters["fluxid"]
        transmissionid = parameters["transmissionid"]
        xiaimage = xiaedf.xiaimage_number(path, radix, number)
        converter = xiaedftonexus.Converter(
            nxentry=self.temp_nxentry,
            qxrfgeometry=self.qxrfgeometry,
            include_counters=include_counters,
            exclude_counters=exclude_counters,
            instrument=instrument,
            fluxid=fluxid,
            transmissionid=transmissionid,
        )
        nxentry = converter(xiaimage)

    @property
    def name(self):
        parameters = self.parameters
        radix, number = parameters["radix"], parameters["number"]
        return radix + ".{}".format(number)

    @staticmethod
    def _rematch_func(redict):
        method = redict.get("method", "regex")
        if method == "regex":
            return lambda ctrname: re.match(redict["pattern"], ctrname)
        elif method == "equal":
            return lambda ctrname: redict["value"] == ctrname
        else:
            return lambda ctrname: False

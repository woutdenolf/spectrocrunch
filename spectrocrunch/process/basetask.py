# -*- coding: utf-8 -*-

from copy import deepcopy
from abc import ABCMeta, abstractmethod
from future.utils import with_metaclass
import logging
from contextlib import contextmanager
import traceback

from ..io import nxfs
from ..io import target
from ..utils import instance
from ..utils import timing
from ..utils.signalhandling import HandleTermination

logger = logging.getLogger(__name__)


class TaskException(Exception):
    pass


class MissingParameter(TaskException):
    pass


class ParameterError(TaskException):
    pass


class Task(with_metaclass(ABCMeta, object)):
    """Task who's output is a single Nexus HDF5 path"""

    def __init__(self, dependencies=None, outputparent=None, **parameters):
        """
        Args:
            dependencies(Optional(Task or h5fs.Path or str))
            outputparent(Optional(h5fs.Path or str))
        """
        self.outputparent = outputparent
        self.dependencies = dependencies
        self.required_parameters = set()
        self.optional_parameters = set()
        self.parameters = parameters
        if not self.hasdependencies and self.outputparent is None:
            raise ValueError(
                'Specify "outputparent" when task does not have dependencies'
            )

    def __str__(self):
        return "Task '{}'".format(self.output)

    def __repr__(self):
        return self.__str__()

    def run(self):
        """Creates the output atomically"""
        if self.dependencies_done:
            self._ensure_outputdeviceparent()
            if self.exists:
                logger.info("{} already done".format(self))
            else:
                logger.info("{} started ...".format(self))
                with timing.timeit_logger(logger, name=str(self)):
                    with self._atomic_context():
                        self._execute()
        else:
            logger.warning("{} not executed (missing dependencies)".format(self))

    def _atomic_context(self):
        return HandleTermination(
            setup=self._atomic_context_enter, teardown=self._atomic_context_exit
        )

    def _atomic_context_enter(self):
        pass

    def _atomic_context_exit(self, exc_type, exc_value, exc_traceback):
        if exc_type:
            logger.error(
                "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
            )
        return 1  # Exception is handled (do not raise it)

    @property
    def output(self):
        """
        returns: h5fs.Path
        """
        return self.outputparent[self.outputname]

    @property
    def dependencies_done(self):
        for dependency in self.dependencies:
            if hasattr(dependency, "output"):
                if not dependency.done:
                    return False
            else:
                if not dependency.exists:
                    return False
        return True

    @property
    def ready_to_run(self):
        return self.dependencies_done

    @property
    def exists(self):
        return self.output.exists

    @property
    def done(self):
        return self.exists

    @property
    def checksum(self):
        return target.calc_checksum(self.dependencies, self.parameters)

    @property
    def parameters(self):
        return self._parameters

    @parameters.setter
    def parameters(self, value):
        self._parameters = deepcopy(value)
        # Default parameters + mark required/optional
        self._parameters_defaults()
        # Missing required parameters?
        for p in self.required_parameters:
            if p not in self._parameters:
                raise MissingParameter("{} of {}".format(p, type(self)))
        # Remove unknown parameters
        # Dev: subclass could overwrite _parameters_filter
        #      and return None to skip this step
        lst = self._parameters_filter()
        if lst:
            self._parameters = {k: v for k, v in self._parameters.items() if k in lst}

    @property
    def dependencies(self):
        """
        Returns:
            list(Task or h5fs.Path)
        """
        return self._dependencies

    @dependencies.setter
    def dependencies(self, value):
        """
        Args:
            value(list or Task or h5fs.Path or str)
        """
        if instance.isstringarray(value):
            self._dependencies = [nxfs.factory(v) for v in value]
        elif instance.isarray(value):
            self._dependencies = value
        else:
            if value is None:
                self._dependencies = []
            elif instance.isstring(value):
                self._dependencies = [nxfs.factory(value)]
            else:
                self._dependencies = [value]
        self.default_outputparent = None

    @property
    def default_outputparent(self):
        if self._default_outputparent is None:
            try:
                previous = self.previous_outputs
            except AttributeError:
                pass
            else:
                if previous:
                    self._default_outputparent = previous[-1].parent
                else:
                    self._default_outputparent = None
        return self._default_outputparent

    @default_outputparent.setter
    def default_outputparent(self, value):
        self._default_outputparent = value

    @property
    def previous_outputs(self):
        """
        Returns:
            list(h5fs.Path)
        """
        return [out for out in self.previous_outputs_iter()]

    def previous_outputs_iter(self):
        """
        Yields h5fs.Path
        """
        for dependency in self.dependencies:
            if isinstance(dependency, Task):
                yield dependency.output
            else:
                yield dependency

    @property
    def hasdependencies(self):
        return bool([x for x in self.dependencies])

    def _parameters_defaults(self):
        self.parameters["name"] = self.parameters.get("name", self.method)
        self.required_parameters |= {"name", "method"}

    def _parameters_filter(self):
        return self.required_parameters | self.optional_parameters

    @property
    def method(self):
        return self.parameters["method"]

    @property
    def outputname(self):
        return self.parameters["name"]

    @property
    def outputparent(self):
        if self._outputparent is None:
            return self.default_outputparent
        else:
            return self._outputparent

    @outputparent.setter
    def outputparent(self, value):
        if instance.isstring(value):
            value = nxfs.Path(value)
        self._outputparent = value

    def _ensure_outputdeviceparent(self):
        device = self.outputparent.device
        if device:
            device.parent.mkdir()
        else:
            self.outputparent.mkdir()

    def find_dependency(self, **parameters):
        """
        Returns:
            nxfs._NXprocess or None
        """
        for nxprocess in self.previous_outputs_iter():
            if nxprocess.is_nxclass(u"NXprocess"):
                config = nxprocess.config.read()
                if all(config.get(k, None) == v for k, v in parameters.items()):
                    return nxprocess
        return None

    @abstractmethod
    def _execute(self):
        pass

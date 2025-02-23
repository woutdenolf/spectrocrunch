import h5py
import numpy as np
import contextlib
from datetime import datetime
import json
import logging
import dateutil.tz
import dateutil.parser
import traceback
from importlib.metadata import version as get_version

from . import fs
from . import h5fs
from ..utils import units
from ..utils import instance
from ..utils import hashing
from ..patch import jsonpickle
from . import target

PROGRAM_NAME = "spectrocrunch"
PROGRAM_VERSION = get_version("spectrocrunch")
DEFAULT_PLOT_NAME = "defaultplot"

logger = logging.getLogger(__name__)


class NexusException(h5fs.FileSystemException):
    """
    Base class for generic Nexus exceptions.
    """

    pass


class NexusFormatException(NexusException):
    """
    Raised when the hdf5 structure doesn't match the Nexus standards.
    """

    pass


class NexusProcessWrongHash(NexusException):
    """
    Raised when the NXprocess exists with a different hash
    """

    pass


tzlocal = dateutil.tz.tzlocal()


def now():
    return datetime.now(tzlocal)


def timestamp():
    return now().isoformat()


def parse_timestamp(tm):
    try:
        return dateutil.parser.parse(tm)
    except ValueError:
        pass
    return tm


class Path(h5fs.Path):
    def mkfile(self, data=None, units=None, attributes=None, **params):
        if attributes is None:
            attributes = {}
        if instance.isquantity(data):
            if units is not None:
                data = data.to(units)
            units = data.units
            data = data.magnitude
        if units is not None:
            if not instance.isstring(units):
                units = "{:~}".format(units)
            attributes["units"] = units
        return super(Path, self).mkfile(data=data, attributes=attributes, **params)

    @property
    def factorycls(self):
        return Path

    def factory(self, path):
        if not isinstance(path, h5fs.Path):
            path = super(Path, self).factory(path)
        cls = None
        nxclass = path.nxclass
        if nxclass == "NXdata":
            cls = _NXdata
        elif nxclass == "NXprocess":
            cls = _NXprocess
        elif nxclass == "NXnote":
            cls = _NXnote
        elif nxclass == "NXmonochromator":
            cls = _NXmonochromator
        elif nxclass == "NXcollection":
            if path.name == "positioners" and path.parent.nxclass == "NXinstrument":
                cls = _Positioners
            else:
                cls = _NXcollection
        if cls is None:
            return path
        else:
            return cls(path, h5file=self._h5file)

    @property
    def nxclass(self):
        if self.exists:
            try:
                return instance.asunicode(self.get_stat("NX_class", default=None))
            except Exception:
                logger.warning("Corrupted " + str(self))
                logger.warning(traceback.format_exc())
        return None

    def updated(self):
        with self.h5open():
            tm = timestamp()
            for path in self.iterup_is_nxclass("NXentry", "NXsubentry"):
                path["end_time"].write(mode="w", data=tm)
            for path in self.iterup_is_nxclass("NXprocess", "NXnote"):
                path["date"].write(mode="w", data=tm)
            self.nxroot().update_stats(file_update_time=tm)

    def _read_time(self, name, *nxclasses):
        path = self.findfirstup_is_nxclass(*nxclasses)
        if path is not None:
            if name in path:
                return parse_timestamp(path[name].read())
        return None

    @property
    def end_time(self):
        return self._read_time("end_time", "NXentry", "NXsubentry")

    @property
    def start_time(self):
        return self._read_time("start_time", "NXentry", "NXsubentry")

    @property
    def date(self):
        return self._read_time("date", "NXprocess", "NXnote")

    def _init_nxclass(
        self, path, nxclass, nxattributes=None, nxfiles=None, **openparams
    ):
        with self.h5open(**openparams):
            if not path.exists:
                path = path.mkdir()
            nxclassattr = path.get_stat("NX_class", None)
            if nxclassattr:
                if nxclassattr != nxclass:
                    raise fs.AlreadyExists(
                        "{} with the wrong class ({} instead of {})".format(
                            path, nxclassattr, nxclass
                        )
                    )
            else:
                if nxattributes:
                    if instance.iscallable(nxattributes):
                        nxattributes = nxattributes()
                    path.update_stats(**nxattributes)
                path.update_stats(NX_class=nxclass)
                if nxfiles:
                    if instance.iscallable(nxfiles):
                        nxfiles = nxfiles()
                    for name, data in nxfiles.items():
                        path[name].mkfile(data=data, **openparams)
                path.updated()
            return self.factory(path)

    def is_nxclass(self, *nxclasses):
        nxclass = self.nxclass
        if nxclasses:
            return any(nxclass == instance.asunicode(cls) for cls in nxclasses)
        else:
            return bool(nxclass)

    def is_not_nxclass(self, *nxclasses):
        nxclass = self.nxclass
        if nxclasses:
            return all(nxclass != instance.asunicode(cls) for cls in nxclasses)
        else:
            return not bool(nxclass)

    def _raise_if_class(self, *nxclasses):
        if self.is_nxclass(*nxclasses):
            raise NexusFormatException(
                "NX_class shouldn't be one of these: {}".format(nxclasses)
            )

    def _raise_ifnot_class(self, *nxclasses):
        if self.is_not_nxclass(*nxclasses):
            raise NexusFormatException(
                "NX_class should be one of these: {}".format(nxclasses)
            )

    def iterup_is_nxclass(self, *nxclasses):
        with self.h5open():
            for path in self.iterup:
                if path.is_nxclass(*nxclasses):
                    yield self.factory(path)

    def iterup_isnot_nxclass(self, *nxclasses):
        with self.h5open():
            for path in self.iterup:
                if path.is_not_nxclass(*nxclasses):
                    yield self.factory(path)

    def iter_is_nxclass(self, *nxclasses):
        with self.h5open():
            for path in self:
                if path.is_nxclass(*nxclasses):
                    yield self.factory(path)

    def iter_nxentry(self):
        for entry in self.root.iter_is_nxclass("NXentry"):
            yield entry

    def iter_nxprocess(self, searchallentries=False):
        if searchallentries:
            entry = None
        else:
            entry = self.findfirstup_is_nxclass("NXentry")
        if entry is None:
            for entry in self.iter_nxentry():
                for process in entry.iter_is_nxclass("NXprocess"):
                    yield process
        else:
            for process in entry.iter_is_nxclass("NXprocess"):
                yield process

    def iter_application(self, searchallentries=False):
        if searchallentries:
            entry = None
        else:
            entry = self.findfirstup_is_nxclass("NXentry")
        if entry is None:
            for entry in self.iter_nxentry():
                for subentry in entry.iter_is_nxclass("NXsubentry"):
                    if "definition" in subentry:
                        yield subentry
        else:
            for subentry in entry.iter_is_nxclass("NXsubentry"):
                if "definition" in subentry:
                    yield subentry

    def findfirstup_is_nxclass(self, *nxclasses):
        path = None
        for path in self.iterup_is_nxclass(*nxclasses):
            return path
        return None

    def findfirstup_is_name(self, *names):
        path = None
        for path in self.iterup:
            if path.name in names:
                return path
        return None

    def nxroot(self, **openparams):
        return self._init_nxclass(
            self.root, "NXroot", nxattributes=self._nxattributes_nxroot, **openparams
        )

    def _nxattributes_nxroot(self):
        return {
            "file_name": self._h5file.path,
            "file_time": timestamp(),
            "HDF5_Version": h5py.version.hdf5_version,
            "h5py_version": h5py.version.version,
            "creator": PROGRAM_NAME,
        }

    def nxentry(self, name=None, **openparams):
        """Get NXentry (first looks in parents) and create if it does not exist"""
        path = self.findfirstup_is_nxclass("NXentry")
        if path is not None:
            if path.name == name or name is None:
                return path

        root = self.nxroot(**openparams)
        if name is None:
            for path in self.iterup:
                if path.parent == root:
                    break
            if path == root:
                name = self.next_nxentry_name()
            else:
                name = path.name

        return self._init_nxclass(
            root[name], "NXentry", nxfiles=self._nxfiles_nxentry, **openparams
        )

    def last_nxentry(self):
        entry = None
        for entry in self.root.iter_is_nxclass("NXentry"):
            pass
        return entry

    def new_nxentry(self, entry=None):
        name = self.next_nxentry_name(entry=entry)
        return self.nxentry(name=name)

    def next_nxentry_name(self, entry=None):
        """New entry name based on the naming logic of a given entry

        Args:
            entry(Optional(Path)): last entry when missing

        Returns:
            str
        """
        # If argument not an NXentry: find it in the parents
        if entry is not None:
            if entry.is_not_nxclass("NXentry"):
                entry = entry.findfirstup_is_nxclass("NXentry")

        # Default NXentry
        if entry is None:
            entry = self.last_nxentry()

        # Next entry name
        if entry is None:
            name = "entry0001"
        else:
            name = target.Name(entry.name)
            parent = entry.parent
            while parent[str(name)].exists:
                name += 1
            name = str(name)
        return name

    def nxsubentry(self, name, **openparams):
        self._raise_ifnot_class("NXentry", "NXsubentry")
        return self._init_nxclass(
            self[name], "NXsubentry", nxfiles=self._nxfiles_nxentry, **openparams
        )

    def _nxfiles_nxentry(self):
        return {"start_time": timestamp()}

    def nxdata(self, name, **openparams):
        self._raise_if_class("NXroot")
        return self._init_nxclass(self[name], "NXdata", **openparams)

    def nxnote(self, name, **openparams):
        self._raise_if_class("NXroot")
        return self._init_nxclass(
            self[name], "NXnote", nxfiles=self._nxfiles_nxnote, **openparams
        )

    def _nxfiles_nxnote(self):
        return {"date": timestamp()}

    def nxprocess(
        self,
        name,
        parameters=None,
        dependencies=None,
        searchallentries=True,
        noincrement=False,
        **openparams,
    ):
        """Creates NXentry and NXprocess when needed.

        Args:
            name(str): only used when creating a new process
            parameters(dict):
            dependencies(list): a list of object with the "checksum" attribute
            searchallentries(bool)
            noincrement(bool)
            openparams(Optional(dict)): overwrite device open parameters (if allowed)

        Returns:
            Path: exists and matching name (not necessarily equal)
        """
        # Existing process with matching name
        process = self.find_nxprocess(
            name=name,
            parameters=parameters,
            dependencies=dependencies,
            searchallentries=searchallentries,
        )
        if process is not None:
            return process

        # Get/create NXentry
        entry = self.nxentry()

        # New process path
        path = entry[name]
        if path.exists:
            if noincrement or not path.is_nxclass("NXprocess"):
                raise fs.AlreadyExists(path)
            name = target.Name(name)
            path = entry[str(name)]
            while path.exists:
                name += 1
                path = entry[str(name)]

        # Atomically create NXprocess with its parameters
        process = self._init_nxclass(
            path, "NXprocess", nxfiles=self._nxfiles_nxprocess, **openparams
        )
        process.set_config(parameters, dependencies=dependencies)
        return process

    def find_nxprocess(
        self, name=None, parameters=None, dependencies=None, searchallentries=False
    ):
        """
        Get the path of an NXprocess, identified by its parameters
        and dependencies (and optionally the name).

        Args:
            name(str):
            parameters(dict):
            dependencies(list): a list of object with the "checksum" attribute

        Returns:
            Path or None
        """
        if name:
            match = target.match(name)
        else:

            def match(x):
                return True

        checksum = target.calc_checksum(dependencies, parameters)
        for process in self.iter_nxprocess(searchallentries=searchallentries):
            if match(process.name):
                if process.checksum == checksum:
                    return process
        return None

    def find_application(self, entryname=None, definition=None):
        for application in self.iter_application(searchallentries=True):
            if entryname:
                if application.parent.name != entryname:
                    continue
            if definition:
                if application["definition"].read() != definition:
                    continue
            return application
        return None

    def _nxfiles_nxprocess(self):
        return {
            "program": PROGRAM_NAME,
            "version": PROGRAM_VERSION,
            "date": timestamp(),
        }

    def nxinstrument(self, **openparams):
        path = self.findfirstup_is_nxclass("NXinstrument")
        if path is not None:
            return path
        entry = self.nxentry()
        return self._init_nxclass(entry["instrument"], "NXinstrument", **openparams)

    def measurement(self, nxclass="NXcollection", **openparams):
        path = self.findfirstup_is_name("measurement")
        if path is not None:
            if path.nxclass:
                return path
        entry = self.nxentry()
        return self._init_nxclass(entry["measurement"], nxclass, **openparams)

    def nxcollection(self, name, nxattributes=None, nxfiles=None, **openparams):
        self._raise_if_class("NXroot")
        return self._init_nxclass(
            self[name],
            "NXcollection",
            nxattributes=nxattributes,
            nxfiles=nxfiles,
            **openparams,
        )

    def nxdetector(self, name, nxattributes=None, nxfiles=None, **openparams):
        instrument = self.nxinstrument(**openparams)
        return self._init_nxclass(
            instrument[name],
            "NXdetector",
            nxattributes=nxattributes,
            nxfiles=nxfiles,
            **openparams,
        )

    def nxmonochromator(
        self, name="monochromator", nxattributes=None, nxfiles=None, **openparams
    ):
        instrument = self.nxinstrument(**openparams)
        return self._init_nxclass(
            instrument[name],
            "NXmonochromator",
            nxattributes=nxattributes,
            nxfiles=nxfiles,
            **openparams,
        )

    def application(
        self, name, definition=None, nxattributes=None, nxfiles=None, **openparams
    ):
        _ = self.nxentry(**openparams)
        if definition is None:
            definition = name
        if not nxfiles:
            nxfiles = {}
        nxfiles.update(self._nxfiles_nxentry())
        nxfiles["definition"] = nxfiles.get("definition", definition)
        return self._init_nxclass(
            self[name],
            "NXsubentry",
            nxattributes=nxattributes,
            nxfiles=nxfiles,
            **openparams,
        )

    def positioners(self, **openparams):
        path = self.findfirstup_is_nxclass("NXprocess")
        if path is None:
            path = self.nxinstrument(**openparams)
        else:
            path = path.results
        return path.nxcollection("positioners", **openparams)

    def mark_default(self):
        path = self
        nxdata = None
        for parent in self.parent.iterup:
            nxclass = path.nxclass
            parentnxclass = parent.nxclass
            if parentnxclass == "NXdata":
                parent.default_signal(path.name)
            elif nxclass == "NXentry":
                parent.update_stats(default=path.name)
            elif parentnxclass is not None:
                if nxclass == "NXdata":
                    nxdata = path
                if nxdata:
                    # if parentnxclass == u'NXentry' and nxclass != u'NXdata':
                    #    default = parent[DEFAULT_PLOT_NAME]
                    #    if default.islink:
                    #        default.remove()
                    #    dest = parent.relpath(nxdata.path)
                    #    nxdata = default.link(dest)
                    dest = parent.relpath(nxdata.path)
                    parent.update_stats(default=dest)
                if parentnxclass == "NXroot":
                    break
            path = parent

    @property
    def default(self):
        path = self.root
        default = path.get_stat("default", default=None)
        while default:
            path = path[default]
            default = path.get_stat("default", default=None)
        return path

    @property
    def dependencies(self):
        self._raise_ifnot_class("NXentry", "NXprocess")
        deproot = self.nxcollection("dependencies")
        order = deproot.get_stat("order", [])
        return [deproot[name] for name in order]

    def add_dependency(self, dependency):
        self._raise_ifnot_class("NXentry", "NXprocess")
        if instance.isstring(dependency):
            dependency = factory(dependency)
        deproot = self.nxcollection("dependencies")
        order = deproot.get_stat("order", None)
        if order is None:
            order = []
        else:
            order = order.tolist()
        name = deproot.nonexisting_name(dependency.name)
        dependency = deproot[name].link(dependency)
        order.append(name)
        deproot.update_stats(order=order)


class _NXPath(Path):
    @contextlib.contextmanager
    def _verify(self):
        with self.h5open():
            self._verify_func()
            yield

    def _verify_func(self):
        return self._raise_ifnot_class(self.NX_CLASS)


class _NXprocess(_NXPath):

    NX_CLASS = "NXprocess"

    @property
    def config(self):
        with self._verify():
            return self.nxnote("configuration")

    @property
    def configpath(self):
        return self["configuration"]

    @property
    def results(self):
        with self._verify():
            return self.nxcollection("results")

    @property
    def resultspath(self):
        return self["results"]

    @property
    def plotselect(self):
        with self._verify():
            return self.nxdata(DEFAULT_PLOT_NAME)

    def set_config(self, parameters, dependencies=None):
        with self._verify():
            # Save parameters
            self.config.write(parameters, typ="jsonpickle")

            # Links to dependencies
            if dependencies:
                for dependency in dependencies:
                    self.add_dependency(dependency)

            # Other info
            self["sequence_index"].write(data=self.sequence_index)
            self.configpath.update_stats(checksum=self.checksum)
            self.updated()

    @property
    def checksum(self):
        with self._verify():
            checksum = self.configpath.get_stat("checksum")
            if checksum is None:
                checksum = target.calc_checksum(
                    self.dependencies, self.confighash, alreadyhash=True
                )
            else:
                checksum = checksum.encode("ascii")

            return checksum

    def valid_checksum(self):
        checksum1 = self.checksum
        checksum2 = target.calc_checksum(
            self.dependencies, self.confighash, alreadyhash=True
        )
        return checksum1 == checksum2

    @property
    def confighash(self):
        if self.configpath.exists:
            return hashing.calchash(self.config.read(parse=True))
        else:
            return None

    @property
    def sequence_index(self):
        """Indexing starts from 1"""
        path = self["sequence_index"]
        if path.exists:
            return path.read()
        else:
            ind = 0
            with self._verify():
                for dependency in self.dependencies:
                    try:
                        ind = max(ind, dependency.sequence_index)
                    except AttributeError:
                        pass
            return ind + 1


class _NXnote(_NXPath):

    NX_CLASS = "NXnote"

    @property
    def mimetype(self):
        if self["type"].exists:
            return self["type"].read()
        else:
            return "text/plain"

    def read(self, parse=True):
        with self._verify():
            data = self["data"].read()
            if parse:
                if self.mimetype == "application/json":
                    data = json.loads(data)
                elif self.mimetype == "application/jsonpickle":
                    data = jsonpickle.loads(data)
            return data

    def _write(self, data, mimetype):
        with self._verify():
            self["type"].mkfile(data=mimetype)
            self["data"].mkfile(data=data)

    def write(self, data, typ="jsonpickle"):
        if instance.isstring(data):
            self._write(data, "text/plain;charset=utf-8")
        else:
            if typ == "json":
                self._write(json.dumps(data, indent=2), "application/json")
            elif typ == "jsonpickle":
                self._write(jsonpickle.dumps(data, indent=2), "application/jsonpickle")
            else:
                self._write(str(data), "text/plain;charset=utf-8")


class _NXdata(_NXPath):

    NX_CLASS = "NXdata"

    @property
    def signal(self):
        with self._verify():
            with self.open():
                name = self.get_stat("signal", None)
                if name:
                    return self[name]
                else:
                    return None

    def _signal_props(self, prop, default):
        signal = self.signal
        if signal:
            with signal.open(mode="r") as dset:
                if prop == "shape":
                    return dset.shape
                elif prop == "dtype":
                    return dset.dtype
                elif prop == "ndim":
                    return dset.ndim
                elif prop == "size":
                    return dset.size
                elif prop == "interpretation":
                    interpretation = signal.get_stat("interpretation", None)
                    if not interpretation:
                        interpretation = self._interpretation(dset)
                    return interpretation
        return default

    @property
    def shape(self):
        return self._signal_props("shape", ())

    @property
    def dtype(self):
        return self._signal_props("dtype", None)

    @property
    def ndim(self):
        return self._signal_props("ndim", None)

    @property
    def size(self):
        return self._signal_props("size", None)

    @property
    def interpretation(self):
        return self._signal_props("interpretation", None)

    def stackdim(self, default=None):
        interpretation = self.interpretation
        if interpretation == "spectrum":
            return self.ndim - 1
        elif interpretation == "image":
            return 0
        else:
            return default

    @property
    def signals(self):
        with self._verify():
            yield self.signal
            for signal in self.auxiliary_signals:
                yield signal

    @property
    def auxiliary_signals(self):
        with self._verify():
            names = self.get_stat("auxiliary_signals", [])
            for name in names:
                yield self[name]

    def add_signal(
        self,
        name=None,
        path=None,
        title=None,
        interpretation=None,
        units=None,
        **createparams,
    ):
        """
        args:
            name(str)
            signal(Path)
            title(str)
            createparams(dict)
        """
        with self._verify():
            # Check dimensions
            ashape = self.axes_shape
            if not ashape:
                ashape = self.shape
            if ashape:
                if "shape" in createparams:
                    dshape = tuple(createparams["shape"])
                elif "data" in createparams:
                    dshape = createparams["data"].shape
                else:
                    if path is None:
                        spath = self[name]
                    else:
                        spath = path
                    with spath.open(mode="r") as dset:
                        dshape = dset.shape
                if dshape != ashape:
                    raise NexusFormatException(
                        "Data dimensions {} do not match axes dimensions {}".format(
                            dshape, ashape
                        )
                    )

            # Create the signal dataset
            if path:
                if name is None:
                    name = path.name
                self[name].link(path)
            elif name:
                if title is None:
                    title = name
                signalpath = self[name]
                with signalpath.open(**createparams) as node:
                    if interpretation is None:
                        interpretation = self._interpretation(node)
                    attrs = {"long_name": title}
                    if interpretation:
                        attrs["interpretation"] = interpretation
                    if units is not None:
                        attrs["units"] = units
                    signalpath.update_stats(**attrs)
            else:
                raise ValueError('Provide either "name" or "signal"')

            # Add signal name to NXdata attributes
            with self.open() as node:
                signals = []
                add = self.get_stat("auxiliary_signals", None)
                if add is not None:
                    signals += add.tolist()
                add = self.get_stat("signal", None)
                if add is not None:
                    signals.append(add)
                if signals:
                    self.update_stats(auxiliary_signals=signals)
                self.update_stats(signal=name)

            self.updated()
            return self[name]

    def _interpretation(self, node):
        ndim = node.ndim
        if ndim == 0:
            return "scalar"
        elif ndim == 1:
            return "spectrum"
        else:
            return "image"

    def default_signal(self, name):
        with self._verify():
            with self.open():
                if name == self.get_stat("signal"):
                    return True
                signals = []
                add = self.get_stat("signal", None)
                if add is not None:
                    signals.append(add)
                add = self.get_stat("auxiliary_signals", None)
                if add is not None:
                    signals += add.tolist()
                if name not in signals:
                    return False
                aux = [signal for signal in signals if signal != name]
                if aux:
                    self.update_stats(auxiliary_signals=aux)
                self.update_stats(signal=name)
            self.updated()
        return True

    def remove_signal(self, name):
        with self._verify():
            self[name].remove()
            with self.open():
                signals = self.pop_stat("auxiliary_signals")
                if signals is not None:
                    if self.signal.name == name:
                        aux = signals[:-1]
                        if aux.size:
                            self.update_stats(auxiliary_signals=aux)
                        self.update_stats(signal=signals[-1])
                    else:
                        aux = [signal for signal in signals if signal != name]
                        if aux:
                            self.update_stats(auxiliary_signals=aux)
            self.updated()

    def sort_signals(self, other=None):
        with self._verify():
            with self.open():
                signals = sorted(sig.name for sig in self.signals)
                if other:
                    other_signals = list(sig.name for sig in other.signals)
                    if sorted(other_signals) != signals:
                        return
                    signals = other_signals
                aux = signals[1:]
                if aux:
                    self.update_stats(auxiliary_signals=aux)
                self.update_stats(signal=signals[0])

            self.updated()

    def set_axes(self, *axes):
        """
        Args:
            *axes(3-tuple|str): (name(str),value(array|Path|None),attr(dict|None)) or name(str)
        """
        with self._verify():
            with self.open():
                axes, positioners = self._axes_parse(axes)

                # Check dimensions
                dshape = self.shape
                if dshape:
                    ashape = tuple(
                        self._axis_size(name, value, positioners)
                        for name, value, _ in axes
                    )
                    if dshape != ashape:
                        raise NexusFormatException(
                            "Axes dimensions {} do not match signal dimensions {}".format(
                                ashape, dshape
                            )
                        )

                # Create axes if needed
                names = []
                for name, value, attr in axes:
                    axis = positioners.add_axis(name, value, **attr)
                    if not self[name].exists:
                        self[name].link(axis, soft=True)
                    names.append(name)
                self.update_stats(axes=names)
                self.updated()
                return axes

    def _axes_parse(self, axesin):
        positioners = self.positioners()

        # Axes as list of tuples
        axes = []
        for axis in axesin:
            if instance.isstring(axis):
                name, value, attr = axis, None, {}
            else:
                name, value, attr = axis
                if attr is None:
                    attr = {}
            axes.append((name, value, attr))

        return axes, positioners

    def _axis_size(self, name, value, positioners):
        if isinstance(value, h5fs.Path):
            with value.open(mode="r") as anode:
                return len(anode)
        try:
            return len(value)
        except TypeError:
            return positioners.axis_size(name)

    @property
    def axes(self):
        with self._verify():
            with self.open():
                ret = []
                for name in self.get_stat("axes", []):
                    axispath = self[name]
                    with axispath.open(mode="r") as axis:
                        attrs = {
                            "title": axispath.get_stat("long_name", None),
                            "units": axispath.get_stat("units", None),
                        }
                        values = axis[()]
                        values = units.Quantity(values, attrs["units"])
                        ret.append((name, values, attrs))
                return ret

    @property
    def axes_shape(self):
        with self._verify():
            with self.open():
                ret = []
                for name in self.get_stat("axes", []):
                    with self[name].open(mode="r") as axis:
                        ret.append(axis.size)
                return tuple(ret)


class _NXcollection(_NXPath):

    NX_CLASS = "NXcollection"

    def add_axis(self, name, value, title=None, units=None):
        with self._verify():
            axis = self[name]
            if axis.exists:
                if value is not None:
                    if not self._axis_equal(axis, value):
                        raise ValueError("Axis {} already exists".format(name))
            else:
                if value is None:
                    raise ValueError(
                        "Axis {} does not exist yet and needs data".format(name)
                    )
                if isinstance(value, h5fs.Path):
                    axis.link(value)
                else:
                    if instance.isquantity(value):
                        units = value.units
                        value = value.magnitude
                    if units:
                        try:
                            units = "{:~}".format(units)
                        except ValueError:
                            pass
                    if title is None:
                        if units:
                            title = "{} ({})".format(name, units)
                        else:
                            title = name
                    axis.mkfile(data=value)
                    axis.update_stats(long_name=title)
                    if units:
                        axis.update_stats(units=units)
                self.updated()
            return axis

    def get_axis(self, name):
        with self._verify():
            axis = self[name]
            if axis.exists:
                with axis.open(mode="r") as node:
                    values = node[()]
                    u = axis.get_stat("units", "dimensionless")
                    values = units.Quantity(values, u)
                return values

    def axis_size(self, name):
        with self._verify():
            axis = self[name]
            if axis.exists:
                with axis.open(mode="r") as node:
                    return node.size
            return None

    def _axis_equal(self, axis1, axis2):
        if isinstance(axis1, h5fs.Path):
            with axis1.open(mode="r") as node1:
                if isinstance(axis2, h5fs.Path):
                    with axis2.open(mode="r") as node2:
                        return np.array_equal(node1, node2)
                else:
                    return np.array_equal(node1, axis2)
        else:
            if isinstance(axis2, h5fs.Path):
                with axis2.open(mode="r") as node2:
                    return np.array_equal(axis1, node2)
            else:
                return np.array_equal(axis1, axis2)


class _Positioners(_NXcollection):

    NX_CLASS = "NXcollection"

    def _verify_func(self):
        self._raise_ifnot_class(self.NX_CLASS)
        if self.name != "positioners":
            NexusFormatException(
                'Name should be "positioners" not "{}"'.format(self.name)
            )


class _NXmonochromator(_NXPath):

    NX_CLASS = "NXmonochromator"

    @property
    def energy(self):
        with self._verify():
            path = self["energy"]
            value = path.read()
            u = path.get_stat("units", "keV")
            return units.Quantity(value, u)

    @energy.setter
    def energy(self, value):
        with self._verify():
            self["energy"].mkfile(data=units.Quantity(value, "keV"))


def factory(path, **params):
    return Path(path, **params).factory(path)

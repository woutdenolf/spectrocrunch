# -*- coding: utf-8 -*-

from ..utils.hashable import CompHashable
from ..utils import instance
from ..io import nxfs
from . import axis


class Group(CompHashable):
    """NXprocess subgroup name wrapper"""

    def __init__(self, groupname):
        name = "counters"
        number = -1
        category = -1

        if instance.isnumber(groupname):
            name = "detector{:02d}".format(groupname)
            number = groupname
            category = 0
        elif instance.isstring(groupname):
            if groupname.startswith("xia"):
                groupname = groupname[3:]
            if groupname.startswith("detector"):
                groupname = groupname[8:]

            if groupname.isdigit():
                number = int(groupname)
                name = "detector{:02d}".format(number)
                category = 0
            elif groupname.startswith("S"):
                number = int(groupname[1:])
                name = "detectorS{:01d}".format(number)
                category = 1
            elif any(s in groupname for s in ["parameters", "concentrations"]):
                name = groupname
                category = 0
            else:
                name = groupname
        elif isinstance(groupname, self.__class__):
            name, number, category = (
                groupname.name,
                groupname.number,
                groupname.category,
            )
        elif groupname:
            raise ValueError("Unexpected detector name {}".format(repr(groupname)))
        self.name = name
        self.number = number
        self.category = category

    @property
    def isdetector(self):
        # e.g. xia00, xiaS0
        return self.category >= 0

    @property
    def issum(self):
        # e.g. xiaS0
        return self.category == 1

    @property
    def issingledetector(self):
        # e.g. xia00
        return self.category == 0

    @property
    def xialabel(self):
        if self.issingledetector:
            return "xia{:02d}".format(self.number)
        elif self.issum:
            return "xiaS{:d}".format(self.number)
        else:
            return None

    @property
    def _repr(self):
        return self.name


def regulargriddata(nxgroup):
    """
    Args:
        nxgroup(nxfs.Path): NXdata or NXprocess

    Returns:
        groups(dict): Group:list(nxfs.Path)
        axes(list(Axis)):
        stackdim(int):
    """
    axes = []
    groups = {}

    if nxgroup.nxclass == "NXdata":
        it = [nxgroup]
        stackdim = None
    elif nxgroup.nxclass == "NXprocess":
        progname = nxgroup["program"].read()
        if progname == nxfs.PROGRAM_NAME:
            stackdim = nxgroup.config.read().get("stackdim", None)
        else:
            raise ValueError(
                'NXprocess from program "{}" is not known'.format(progname)
            )
        it = nxgroup.results.iter_is_nxclass("NXdata")
    else:
        raise ValueError("'{}' should be an NXdata or NXprocess group".format(nxgroup))

    for nxdata in it:
        if nxdata.islink:
            continue
        if stackdim is None:
            stackdim = nxdata.stackdim()
        group = Group(nxdata.name)
        if group in groups:
            raise RuntimeError("Group '{}' appears more than once".format(group))
        axs = [
            axis.factory(values, name=name, title=attrs["title"], type="quantitative")
            for name, values, attrs in nxdata.axes
        ]
        if axes:
            if axes != axs:
                raise RuntimeError("{} has different axes".format(nxdata))
        else:
            axes = axs
        groups[group] = list(nxdata.signals)

    return groups, axes, stackdim

# -*- coding: utf-8 -*-

from ..io import nxfs


def set_default(parent, default):
    it = parent.iter_is_nxclass("NXdata")
    # Default signal
    nxdataselect = None
    for nxdata in it:
        if default is None:
            default = nxdata.signal.name
            nxdataselect = nxdata
        else:
            if nxdata.default_signal(default):
                nxdataselect = nxdata
    # Default NXdata
    if nxdataselect:
        nxdataselect.mark_default()
    else:
        parent.mark_default()
    parent.updated()

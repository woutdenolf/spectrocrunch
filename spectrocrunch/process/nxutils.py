# -*- coding: utf-8 -*-

from ..io import nxfs


def set_default(nxprocess, default):
    it = nxprocess.results.iter_is_nxclass("NXdata")
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
        nxprocess.mark_default()
    nxprocess.updated()

# -*- coding: utf-8 -*-

import logging
from .utils.cli import logging_cliconfig

try:
    from ._version import version as __version__
except ImportError:
    import os

    __version__ = "Local version ({})".format(
        os.path.dirname(os.path.abspath(__file__))
    )

logger = logging.getLogger(__name__)
logging_cliconfig(logger)

# -*- coding: utf-8 -*-

import logging
from ..utils import timing

logger = logging.getLogger(__name__)


def run_sequential(tasks, name=None):
    """
    Args:
        tasks(list(NXtask))
    Returns:
        bool
    """
    with timing.timeit_logger(logger, name=name):
        for task in tasks:
            task.run()
            if not task.output.exists:
                return False
        return True

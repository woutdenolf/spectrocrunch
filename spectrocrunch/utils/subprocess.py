# -*- coding: utf-8 -*-

from __future__ import absolute_import
import os
import errno
import subprocess


def _execute(*args, **kwargs):
    proc = subprocess.Popen(args, **kwargs)
    out, err = proc.communicate()
    return out, err, proc.returncode


def installed(*args):
    try:
        devnull = open(os.devnull)
        _execute(*args, stdout=devnull, stderr=devnull)
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True


def execute(*args, **kwargs):
    try:
        if kwargs.pop('stdout', False):
            kwargs['stdout'] = subprocess.PIPE
        if kwargs.pop('stderr', False):
            kwargs['stderr'] = subprocess.PIPE
        return _execute(*args, **kwargs)
    except OSError as e:
        return None, None, e.errno
# -*- coding: utf-8 -*-

import os
import errno
import shutil
import tempfile
import random
import string


def randomstring(size=6, chars=string.ascii_letters + string.digits):
    # Number of combinations: n^size  (default: 62^6)
    return "".join(random.choice(chars) for _ in range(size))


class Copy(object):
    def __init__(self, filename, copyname):
        self.filename = filename
        self.copyname = copyname

    def __enter__(self):
        shutil.copy2(self.filename, self.copyname)
        return self.copyname

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type:
            if os.path.isfile(self.copyname):
                os.remove(self.copyname)


class TemporaryCopy(object):
    def __init__(self, filename, ext=".tmp"):
        self.filename = filename
        self.tmpfilename = None
        self.ext = ext

    def __enter__(self):
        temp_dir = tempfile.gettempdir()
        temp_name = next(tempfile._get_candidate_names())
        self.tmpfilename = os.path.join(temp_dir, temp_name + self.ext)
        shutil.copy2(self.filename, self.tmpfilename)
        return self.tmpfilename

    def __exit__(self, exc_type, exc_val, exc_tb):
        if os.path.isfile(self.tmpfilename):
            os.remove(self.tmpfilename)
        self.tmpfilename = None


class TemporaryFilename(object):
    def __init__(self, path="", suffix=".tmp", prefix=""):
        if not path:
            path = tempfile.gettempdir()
        self.tmpfilename = temporary_filename(path, suffix=suffix, prefix=prefix)

    def __enter__(self):
        return self.tmpfilename

    def __exit__(self, exc_type, exc_val, exc_tb):
        if os.path.exists(self.tmpfilename):
            os.remove(self.tmpfilename)


def temporary_filename(path, suffix=".tmp", prefix=""):
    if not path:
        path = tempfile.gettempdir()
    return os.path.join(path, prefix + randomstring() + suffix)


def mkdir(path):
    if not path:
        return  # current working directory
    try:
        os.makedirs(path)
    except OSError as e:
        if e.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise e

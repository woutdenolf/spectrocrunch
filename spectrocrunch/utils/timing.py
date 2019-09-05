# -*- coding: utf-8 -*-

import datetime
import sys
import time
import logging
from contextlib import contextmanager


_logger = logging.getLogger(__name__)


def taketimestamp():
    return datetime.datetime.now()


def hms(seconds):
    sec = seconds
    hours = int(sec/3600.)
    sec -= hours*3600
    min = int(sec/60.)
    sec -= min*60
    sec = int(sec)
    return (hours, min, sec)


def strseconds(seconds):
    return "{:d}h {:d}m {:d}s".format(*hms(seconds))


def printtimeelapsed(T0, logger=None, text="Elapsed time"):
    if logger is None:
        logger = _logger
    hours, min, sec = hms((datetime.datetime.now()-T0).seconds)
    logger.info("{}: {:d}h {:d}m {:d}s".format(text, hours, min, sec))


class ProgressLogger(object):
    def __init__(self, logger=None):
        if logger is None:
            logger = _logger
        self.logger = logger
        self.setn(1)
        self.start()

    def start(self):
        self.T0 = datetime.datetime.now()

    def setn(self, nmax):
        self.nmax = nmax
        self.ncurrent = 0
        self.setnfine(1)

    def setnfine(self, nmaxfine):
        self.nmaxfine = nmaxfine
        self.ncurrentfine = 0

    def ndone(self, n):
        self.ncurrent += n
        self.ncurrentfine = 0

    def ndonefine(self, n):
        self.ncurrentfine += n

    def printprogress(self):
        # Ttot = Telapsed + Tleft
        # Telapsed/Ttot = ncurrent/nmax = pdone
        # Tleft = Ttot - Telapsed = Telapsed*(1/pdone - 1)
        nmax = self.nmax*self.nmaxfine
        ncurrent = self.ncurrent*self.nmaxfine + self.ncurrentfine
        Telapsed = (datetime.datetime.now()-self.T0).seconds
        hours, min, sec = hms(Telapsed)
        hours2, min2, sec2 = hms(Telapsed*(nmax/float(ncurrent)-1))
        pdone = int(ncurrent/float(nmax)*100)
        self.logger.info("  Elapsed: %dh %dm %ds (%d%%)  Left: %dh %dm %ds (%d%%)" % (
            hours, min, sec, pdone, hours2, min2, sec2, 100-pdone))
        sys.stdout.flush()


@contextmanager
def timeit(name=""):
    t0 = time.time()
    yield
    t1 = time.time()
    print("Execution time ({}): {}".format(name, t1-t0))


@contextmanager
def timeit_logger(logger=None, name=""):
    T0 = taketimestamp()
    yield
    if name:
        text = "Elapsed time ({})".format(name)
    else:
        text = None
    printtimeelapsed(T0, logger=logger, text=text)

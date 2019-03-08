# -*- coding: utf-8 -*-
#
#   Copyright (C) 2018 European Synchrotron Radiation Facility, Grenoble, France
#
#   Principal author:   Wout De Nolf (wout.de_nolf@esrf.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import argparse
from contextlib import contextmanager
import collections
import logging

from ..utils import instance

logger = logging.getLogger(__name__)


def _asstring(x):
    if instance.isstring(x):
        return "\"{}\"".format(x)
    else:
        return str(x)


def jobname(name, args, more=True):
    if more:
        add = ", ..."
    else:
        add = ""

    return "{}({}{})".format(name, ", ".join([_asstring(a) for a in args]), add)


class Job(object):

    def __init__(self, name, func, args, kwargs):
        if not instance.iscallable(func):
            raise TypeError(
                "Job second argument \"{}\" must be callable".format(func))
        if not instance.isiterable(args):
            raise TypeError(
                "Job third argument \"{}\" must be iterable".format(args))
        if not isinstance(kwargs, dict):
            raise TypeError(
                "Job forth argument \"{}\" must be a dictionary".format(kwargs))
        self.name = name
        self.func = func
        self.args = args
        self.kwargs = kwargs

    def __call__(self):
        logger.info("Execute job: {}".format(self))
        self.func(*self.args, **self.kwargs)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.__str__()


class Jobs(collections.MutableSequence):

    def __init__(self, *args):
        self.list = list()
        self.extend(list(args))

    def __len__(self):
        return len(self.list)

    def __getitem__(self, i):
        return self.list[i]

    def __delitem__(self, i):
        del self.list[i]

    def __setitem__(self, i, v):
        self.list[i] = Job(*v)

    def insert(self, i, v):
        self.list.insert(i, Job(*v))

    def __str__(self):
        fmt = "{{:{}}}. {{}}".format(len(str(max(len(self)-1, 1))))

        return "\n".join([fmt.format(i, item) for i, item in enumerate(self.list, 1)])


@contextmanager
def context():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--jobs", type=int, default=[],
                        nargs='+', help="execute jobs")
    parser.add_argument("-a", "--ajob", type=int,
                        default=0, help="execute jobs from a")
    parser.add_argument("-b", "--bjob", type=int,
                        default=0, help="execute jobs from b")
    parser.add_argument("-n", "--njobs", action="store_true",
                        help="number of jobs")
    parser.add_argument("-l", "--ljobs", action="store_true", help="list jobs")
    args, unknown = parser.parse_known_args()

    jobs = Jobs()

    yield jobs

    if args.njobs:
        print(len(jobs))
    elif args.ljobs:
        print(jobs)
    elif args.jobs or args.ajob > 0 or args.bjob > 0:
        ind = list(args.jobs)
        if args.ajob > 0 or args.bjob > 0:
            if args.ajob > 0:
                a = args.ajob
            else:
                a = 1
            if args.bjob > 0:
                b = args.bjob
            else:
                b = len(jobs)
            a, b = min(a, b), max(a, b)
            ind = ind+range(a, b+1)

        for i in ind:
            try:
                job = jobs[i-1]
            except IndexError:
                logger.error("Invalid job {} of {}".format(i, len(jobs)))
            else:
                job()
    else:
        for job in jobs:
            job()

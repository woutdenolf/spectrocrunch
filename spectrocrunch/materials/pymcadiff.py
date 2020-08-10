# -*- coding: utf-8 -*-

from __future__ import print_function
from pprint import pprint
from deepdiff import DeepDiff
from PyMca5.PyMcaIO import ConfigDict


def diff(cfgfile1, cfgfile2, view="tree", significant_digits=5, **kwargs):
    fconfig1 = ConfigDict.ConfigDict()
    fconfig1.read(cfgfile1)
    fconfig2 = ConfigDict.ConfigDict()
    fconfig2.read(cfgfile2)

    ddiff = DeepDiff(
        fconfig1,
        fconfig2,
        ignore_order=True,
        view=view,
        significant_digits=significant_digits,
        **kwargs
    )
    pprint(ddiff)

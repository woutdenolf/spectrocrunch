# -*- coding: utf-8 -*-

import inspect


def printstack(i):
    frame, filename, line_number, function_name, lines, index = inspect.stack()[i]
    print "Caller: file = {}, function = {}, line {}".format(
        filename, function_name, line_number
    )


def printcaller():
    printstack(3)

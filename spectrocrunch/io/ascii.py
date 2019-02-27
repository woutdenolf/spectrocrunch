# -*- coding: utf-8 -*-
#
#   Copyright (C) 2015 European Synchrotron Radiation Facility, Grenoble, France
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

import os
import re
import numpy as np
from .utils import mkdir


class asciifile(object):
    """Interface to an ascii file
    """

    def __init__(self, filename):
        if not os.path.isfile(filename):
            raise IOError("File %s does not exist." % filename)

        self.filename = filename

    def autoread_alphatable(self, dtype=np.float32):
        """Read an ascii file with following format:
            1. some header
            2. line with column headers
            3. table with numbers
        """
        with open(self.filename, 'r') as f:
            data = f.read()

        # Regular expressions
        number = r'(?:[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)'
        blank = r'[ \t\f]'
        endline = r'[\r\n?|\n]'
        dataline = r'((?:%s*%s)+)%s*(?:%s|\Z)' % (blank,
                                                  number, blank, endline)
        datalines = r'((?:\s*%s)+)\s*' % number
        header = r'(.*?)%s' % endline
        tabel = r'(?:%s|\A)%s%s' % (endline, header, datalines)

        # Find the first occurance of a tabel with header
        # (Repeated captures are not allowed, only returns the last one)
        p = re.search(tabel, data)
        if p is None:
            raise IOError("No table with header found.")
        if p.lastindex != 2:
            raise IOError("No table with header found.")

        # Extract data
        colheader = re.split(r'\s*', p.group(1))
        data = np.array([re.split(r'\s*', line.strip())
                         for line in re.split(endline, p.group(2))], dtype=dtype)
        return (data, colheader)


class Writer(object):

    def __init__(self, filename=None, mode='w'):
        self.filename = filename
        self.mode = mode
        self.open_file = None

    def __enter__(self):
        if self.filename:
            mkdir(os.path.dirname(self.filename))
            self.open_file = open(self.filename, self.mode)
        else:
            self.open_file = None
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.open_file:
            self.open_file.close()
        self.open_file = None

    def write(self, line):
        self.writelines([line])

    def writelines(self, lines):
        if self.open_file:
            self.open_file.write('\n'.join(lines)+'\n')


class Logger(Writer):

    def __init__(self, filename=None, tab=4, **kwargs):
        self.tab = ' '*tab
        super(Logger, self).__init__(filename=filename, **kwargs)

    def writelines(self, lines, tab=0):
        lines = ["{}{}".format(self.tab*tab, line) for line in lines]
        super(Logger, self).writelines(lines)
        for line in lines:
            print(line)

    def write(self, line, tab=0):
        lines = re.split("\r\n|\r|\n", str(line))
        self.writelines(lines, tab=tab)

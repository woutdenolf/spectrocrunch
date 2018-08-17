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

try:
    from ._version import version as __version__
except ImportError:
    import os
    __version__ = "Local version ({})".format(os.path.dirname(os.path.abspath(__file__)))

# Add default logging handlers when non present
import logging
logger = logging.getLogger(__name__)

def logging_cliconfig():
    """Configure logging from command-line options:
        --log=...     Log level
        --logfile=... Log file
        --stdout=...  Redirect stdout to a file
        --stderr=...  Redirect stderr to a file
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--log',default="",type=str,help="Log level")
    parser.add_argument('--logfile',default="",type=str,help="Log file")
    parser.add_argument('--stdout',default="",type=str,help="Log file for what normally goes to stdout")
    parser.add_argument('--stderr',default="",type=str,help="Log file for what normally goes to stderr")
    args = parser.parse_args()
    
    hashandlers = bool(logger.handlers)
    
    if args.log and not hashandlers:
        logger.setLevel(args.log.upper())
    if args.logfile:
        logging_filehandler(args.logfile)
    if args.stdout:
        logging_filehandler(args.stdout,error=False)
    elif not hashandlers:
        logging_stdhandler(error=False)
    if args.stderr:
        logging_filehandler(args.stderr,error=True)
    elif not hashandlers:
        logging_stdhandler(error=True)
        
def logging_addformat(handler):
    #formatter = logging.Formatter("%(asctime)s %(levelname)s [%(module)s]: %(message)s")
    #formatter = logging.Formatter('[%(asctime)s] p%(process)s {%(pathname)s:%(lineno)d} %(levelname)s - %(message)s','%m-%d %H:%M:%S')
    formatter = logging.Formatter("%(levelname)s:%(name)s: %(message)s")
    handler.setFormatter(formatter)

def logging_stdhandler(error=True):
    """Add stdout handler
    """
    import sys
    if error:
        handler = sys.stderr
    else:
        handler = sys.stdout
    handler = logging.StreamHandler(handler)
    logging_configure(handler,error=error)
    logger.addHandler(handler)

def logging_filehandler(filename,error=None):
    """Add file handler
    """
    handler = logging.FileHandler(filename)
    logging_configure(handler,error=error)
    logger.addHandler(handler)

def logging_configure(handler,error=None):
    logging_addformat(handler)
    logging_addfilter(handler,error=error)

def logging_addfilter(handler,error=None):
    if error is None:
        return
        
    if error:
        func = lambda level: level >= logging.WARNING
    else:
        func = lambda level: level < logging.WARNING

    class Filter(logging.Filter):
        def filter(self, record):
            return func(record.levelno)
    
    handler.addFilter(Filter())

logging_cliconfig()
    

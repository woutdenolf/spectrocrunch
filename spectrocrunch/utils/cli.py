"""
Capture CLI options for logging
"""

import logging
import argparse


def logging_cliconfig(logger):
    """Configure logging from command-line options:
    --log=...     Log level
    --logfile=... Log file
    --stdout=...  Redirect stdout to a file
    --stderr=...  Redirect stderr to a file
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--log", default="", type=str, help="Log level")
    parser.add_argument("--logfile", default="", type=str, help="Log file")
    parser.add_argument(
        "--stdout",
        default="",
        type=str,
        help="Log file for what normally goes to stdout",
    )
    parser.add_argument(
        "--stderr",
        default="",
        type=str,
        help="Log file for what normally goes to stderr",
    )
    args, _ = parser.parse_known_args()

    hashandlers = logger_has_handlers(logger)

    if args.log and not hashandlers:
        logger.setLevel(args.log.upper())
    if args.logfile:
        logging_filehandler(args.logfile, logger)
    if args.stdout:
        logging_filehandler(args.stdout, logger, error=False)
    elif not hashandlers:
        logging_stdhandler(logger, error=False)
    if args.stderr:
        logging_filehandler(args.stderr, logger, error=True)
    elif not hashandlers:
        logging_stdhandler(logger, error=True)


def logger_has_handlers(logger):
    while logger is not None:
        if bool(logger.handlers):
            return True
        logger = logger.parent
    return False


def logging_addformat(handler):
    formatter = logging.Formatter("%(levelname)s:%(name)s: %(message)s")
    handler.setFormatter(formatter)


def logging_stdhandler(logger, error=True):
    """Add stdout handler"""
    import sys

    if error:
        handler = sys.stderr
    else:
        handler = sys.stdout
    handler = logging.StreamHandler(handler)
    logging_configure(handler, error=error)
    logger.addHandler(handler)


def logging_filehandler(filename, logger, error=None):
    """Add file handler"""
    handler = logging.FileHandler(filename)
    logging_configure(handler, error=error)
    logger.addHandler(handler)


def logging_configure(handler, error=None):
    logging_addformat(handler)
    logging_addfilter(handler, error=error)


def logging_addfilter(handler, error=None):
    if error is None:
        return
    if error:

        def func(level):
            return level >= logging.WARNING

    else:

        def func(level):
            return level < logging.WARNING

    class Filter(logging.Filter):
        def filter(self, record):
            return func(record.levelno)

    handler.addFilter(Filter())


def getLogger(name, filename):
    if name == "__main__":
        logname = filename
    else:
        logname = name
    logger = logging.getLogger(logname)
    if name == "__main__":
        logging_cliconfig(logger)
    return logger

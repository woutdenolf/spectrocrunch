import signal
import logging
import sys

logger = logging.getLogger(__name__)


def gethandlers(signames):
    ret = {}
    for attr in signames:
        signum = getattr(signal, attr, None)
        try:
            ret[signum] = signal.getsignal(signum)
        except (ValueError, TypeError):
            pass
    return ret


def termination_handlers():
    # SIGTERM: politely ask for termination
    # SIGINT: user asks for termination (CTRL-C)
    # SIGQUIT: SIGINT but do not cleanup
    # SIGKILL: immediate termination
    # SIGHUP: controlling process terminated
    return gethandlers(["SIGTERM", "SIGINT", "SIGHUP"])


def allhandlers():
    return gethandlers([name for name in dir(signal) if name.startswith("SIG")])


def replace_handlers(newhandler, handlers=None):
    if handlers is None:
        handlers = allhandlers()
    oldhandlers = {}
    for signum, oldhander in handlers.items():
        try:
            oldhander = signal.signal(signum, newhandler)
        except (OSError, RuntimeError, ValueError):
            pass
        else:
            oldhandlers[signum] = oldhander
    return oldhandlers


def restore_handlers(oldhandlers):
    for signum, oldhandler in oldhandlers.items():
        signal.signal(signum, oldhandler)


def send_signals(signals, handlers=None):
    if handlers is None:
        handlers = allhandlers()
    for sig, frame in signals:
        handler = handlers.get(sig, None)
        if handler is not None:
            try:
                handler(sig, frame)
            except TypeError:
                pass


def signalname(signum):
    return next(v for v, k in signal.__dict__.items() if k == signum)


class HandleTermination(object):
    """
    Context manager which allows for tear-down on
    exceptions and graceful termination signals
    """

    def __init__(self, setup=None, teardown=None, resend=True):
        self._needcleanup = False
        self._setup = setup
        self._teardown = teardown
        self.resend = resend

    def __enter__(self):
        """Overwrite signal handlers and call user setup"""
        self._needcleanup = True
        self._oldhandlers = replace_handlers(self._newhandler, termination_handlers())
        if self._setup is not None:
            self._setup()

    def _cleanup(self, exc_type, exc_value, exc_traceback):
        """User tear-down and reset signal handlers"""
        ret = 0
        if self._needcleanup:
            if self._teardown is not None:
                ret = self._teardown(exc_type, exc_value, exc_traceback)
            restore_handlers(self._oldhandlers)
            self._needcleanup = False
        return ret

    def _newhandler(self, signum, frame):
        """Converts signal into exception and reset signals.
        Optionally resends the signal after user tear-down.
        """
        msg = "Signal {} received.".format(signalname(signum))
        logger.debug(msg + " Cleanup and resend...")

        # Convert signal to exception for user tear-down function
        try:
            raise RuntimeError(msg)
        except RuntimeError:
            exc_type, exc_value, exc_traceback = sys.exc_info()

        # User tear-down and reset signal handlers
        self._cleanup(exc_type, exc_value, exc_traceback)

        # Resend the signal to the original handlers
        if self.resend:
            send_signals([(signum, frame)], self._oldhandlers)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        """User tear-down and reset signal handlers"""
        return self._cleanup(exc_type, exc_value, exc_traceback)

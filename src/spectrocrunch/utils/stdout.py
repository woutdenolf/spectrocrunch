import sys
import os


class stdout_redirect(object):
    """
    A context manager that redirects stdout for its scope, usage:

    with stdout_redirect():
        os.system('ls -l')
    """

    def __init__(self, to=os.devnull):
        sys.stdout.flush()
        self._origstdout = sys.stdout
        self._oldstdout_fno = os.dup(sys.stdout.fileno())
        self._dest = os.open(to, os.O_WRONLY)

    def __enter__(self):
        self._newstdout = os.dup(1)
        os.dup2(self._dest, 1)
        os.close(self._dest)
        sys.stdout = os.fdopen(self._newstdout, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self._origstdout
        sys.stdout.flush()
        os.dup2(self._oldstdout_fno, 1)

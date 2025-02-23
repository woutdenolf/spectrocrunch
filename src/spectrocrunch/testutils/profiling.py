try:
    import tracemalloc
    import linecache
except ImportError:
    tracemalloc = None
import cProfile
import pstats

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import logging
from contextlib import contextmanager


logger = logging.getLogger(__name__)


def print_malloc_snapshot(snapshot, key_type="lineno", limit=10, units="KB"):
    """
    :param tracemalloc.Snapshot snapshot:
    :param str key_type:
    :param int limit: limit number of lines
    :param str units: B, KB, MB, GB
    """
    n = ["b", "kb", "mb", "gb"].index(units.lower())
    sunits, units = units, 1024**n

    snapshot = snapshot.filter_traces(
        (
            tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
            tracemalloc.Filter(False, "<unknown>"),
        )
    )
    top_stats = snapshot.statistics(key_type)
    total = sum(stat.size for stat in top_stats)

    print("================Memory profile================")
    for index, stat in enumerate(top_stats, 1):
        frame = stat.traceback[0]
        # replace "/path/to/module/file.py" with "module/file.py"
        # filename = os.sep.join(frame.filename.split(os.sep)[-2:])
        filename = frame.filename
        print(
            "#%s: %s:%s: %.1f %s"
            % (index, filename, frame.lineno, stat.size / units, sunits)
        )
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print("    %s" % line)
        if index >= limit:
            break

    other = top_stats[index:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f %s" % (len(other), size / units, sunits))

    print("Total allocated size: %.1f %s" % (total / units, sunits))
    print("============================================")


@contextmanager
def print_malloc_context(**kwargs):
    """
    :param **kwargs: see print_malloc_snapshot
    """
    if tracemalloc is None:
        logger.error("tracemalloc required")
    else:
        tracemalloc.start()
    try:
        yield
    finally:
        if tracemalloc is not None:
            snapshot = tracemalloc.take_snapshot()
            print_malloc_snapshot(snapshot, **kwargs)


@contextmanager
def print_time_context(restrictions=None, sortby="cumtime"):
    pr = cProfile.Profile()
    pr.enable()
    try:
        yield
    finally:
        pr.disable()
        s = StringIO()
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        if restrictions is None:
            restrictions = (0.1,)
        ps.print_stats(*restrictions)
        print("================Time profile================")
        print(s.getvalue())
        print("============================================")


@contextmanager
def profile(memory=True, time=True, memlimit=10, restrictions=None, sortby="cumtime"):
    if not memory and not time:
        yield
    elif memory and time:
        with print_time_context(restrictions=restrictions, sortby=sortby):
            with print_malloc_context(limit=memlimit):
                yield
    elif memory:
        with print_malloc_context(limit=memlimit):
            yield
    else:
        with print_time_context(restrictions=restrictions, sortby=sortby):
            yield

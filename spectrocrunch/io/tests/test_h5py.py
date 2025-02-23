import os
import sys
import errno
import h5py
import logging
import re
import unittest
import multiprocessing
import threading
from testfixtures import TempDirectory

from ...io import excel

logger = logging.getLogger(__name__)


def checkdata(f):
    assert f["data"][()] == 10


def setdata(f):
    f["data"] = 10


def createHolder(Base):
    class Holder(Base):
        def __init__(self, destination, desttype, holdmode, startevent, exitevent):
            super(Holder, self).__init__()
            self.destination = destination
            self.desttype = desttype
            self.holdmode = holdmode
            self.startevent = startevent
            self.exitevent = exitevent

        def run(self):
            try:
                self._run()
            except Exception:
                self._sendevent()
                raise

        def _run(self):
            logger.debug(
                "\nDestination: type = {}, mode = {}".format(
                    self.desttype, self.holdmode
                )
            )

            # Remove destination
            try:
                os.remove(self.destination)
            except OSError as err:
                if err.errno == errno.ENOENT:
                    pass
                # elif err.errno==errno.EISDIR: # does not work on windows
                elif os.path.isdir(self.destination):
                    os.rmdir(self.destination)
                else:
                    print(os.path.exists(self.destination))
                    print(os.path.isdir(self.destination))
                    raise

            # Hold destination
            if self.desttype == "file":
                with open(self.destination, mode="w") as f:
                    f.write("content")
                if self.holdmode:
                    with open(self.destination, mode=self.holdmode) as f:
                        self._sendevent()
                else:
                    self._sendevent()
            elif self.desttype == "dir":
                os.mkdir(self.destination)
                self._sendevent()
            elif self.desttype == "hdf5":
                if self.holdmode:
                    if self.holdmode.startswith("r"):
                        with h5py.File(self.destination, mode="w") as f:
                            setdata(f)
                        with h5py.File(self.destination, mode=self.holdmode) as f:
                            self._sendevent()
                    else:
                        with h5py.File(self.destination, mode=self.holdmode) as f:
                            setdata(f)
                            self._sendevent()
                else:
                    self._sendevent()
            else:
                self._sendevent()

        def _sendevent(self):
            self.startevent.set()
            self.exitevent.wait()

    return Holder


def createClient(Base):
    class Client(Base):
        def __init__(self, destination, mode, output):
            super(Client, self).__init__()
            self.destination = destination
            self.mode = mode
            self.output = output

        def run(self):
            logger.debug("Hdf5 open mode = {}".format(self.mode))
            output = self._testopen(self.mode)
            logger.debug(output)
            updatedata(self.output, output)

        def _testopen(self, mode):
            try:
                with h5py.File(self.destination, mode=mode) as f:
                    output = "OK"
                    if mode == "r":
                        checkdata(f)
            except IOError as err:
                output = errno.errorcode.get(err.errno, None)
                if output is None:
                    errmsg = str(err)
                    m = re.search("errno = ([0-9]+)", errmsg)
                    if m:
                        output = errno.errorcode.get(int(m.groups()[0]), None)
                    if output is None:
                        m = re.search("error message = '(.+)'", errmsg)
                        if m:
                            output = m.groups()[0]
                        elif "file signature not found" in errmsg:
                            output = "H5_SIGERR"
                        elif (
                            "unable to truncate a file which is already open" in errmsg
                        ):
                            output = "H5_TRUNCERR"
                        elif "file exists" in errmsg:
                            output = "H5_EEXIST"
                        elif "file is already open" in errmsg:
                            output = "H5_OPENERR"
                        elif errmsg:
                            output = errmsg
                        else:
                            output = str(err)
            except (KeyError, AssertionError):
                output += "_NODATA"
            return output

    return Client


def datainit(filename):
    try:
        os.remove(filename)
    except OSError as err:
        if err.errno == errno.ENOENT:
            pass


def updatedata(info, value):
    filename = info["filename"]
    sheet_name = info["sheet_name"]
    row = info["row"]
    column = info["column"]
    data = excel.DataFrame.fromexcel(filename, index_col=[0, 1])
    df = data.get(sheet_name, None)
    if df is None:
        df = excel.DataFrame(sheet_name=sheet_name, rowlevels=("Existing", "Lock mode"))
        data[sheet_name] = df
    df.addvalue(row, column, str(value))
    with excel.Writer(filename) as writer:
        for df in data.values():
            df.writer = writer
            df.columnwidth = 25
            df.save()


def run(Base, Event, destination, filename, sheet_name):
    locking = os.environ.get("HDF5_USE_FILE_LOCKING", None)
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

    Holder = createHolder(Base)
    Client = createClient(Base)

    filemode = [None, "r", "r+", "w", "w+", "a", "a+"]
    dirmode = [None]
    lockmode = [None, "r", "r+", "w", "x", "a"]
    nonemode = [None]
    clientmode = ["r", "r+", "w", "x", "a"]

    desttypes = (
        ["file"] * len(filemode)
        + ["dir"] * len(dirmode)
        + ["hdf5"] * len(lockmode)
        + ["None"] * len(nonemode)
    )
    holdmodes = filemode + dirmode + lockmode + nonemode

    startevent = Event()
    exitevent = Event()

    for desttype, holdmode in zip(desttypes, holdmodes):
        for mode in clientmode:
            output = {
                "filename": filename,
                "row": (desttype, str(holdmode)),
                "column": mode,
            }
            holder = Holder(destination, desttype, holdmode, startevent, exitevent)
            holder.start()
            startevent.wait()

            output["sheet_name"] = sheet_name + "_locked"
            client = Client(destination, mode, output)
            client.start()
            client.join()

            exitevent.set()
            holder.join()
            startevent.clear()
            exitevent.clear()

            output["sheet_name"] = sheet_name + "_unlocked"
            client = Client(destination, mode, output)
            client.start()
            client.join()

    if locking:
        os.environ["HDF5_USE_FILE_LOCKING"] = locking


def main(path=None):
    if path is None:
        path = TempDirectory().path
    outfilename = os.path.join(path, "test.xlsx")
    h5filename = os.path.join(path, "test.h5")
    datainit(outfilename)
    run(
        multiprocessing.Process,
        multiprocessing.Event,
        h5filename,
        outfilename,
        "process",
    )
    logger.debug("---------------------------------")
    run(threading.Thread, threading.Event, h5filename, outfilename, "thread")


class test_h5py(unittest.TestCase):
    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()

    @property
    def outfilename(self):
        return os.path.join(self.dir.path, "test.xlsx")

    @property
    def h5filename(self):
        return os.path.join(self.dir.path, "test.h5")

    def test_thread(self):
        datainit(self.outfilename)
        run(
            threading.Thread,
            threading.Event,
            self.h5filename,
            self.outfilename,
            "thread",
        )

    @unittest.skipIf(sys.platform.startswith("win"), "Does not work under Windows")
    def test_process(self):
        datainit(self.outfilename)
        run(
            multiprocessing.Process,
            multiprocessing.Event,
            self.h5filename,
            self.outfilename,
            "process",
        )


def main_test_suite():
    """Test suite including all test suites"""
    testSuite = unittest.TestSuite()
    testSuite.addTest(test_h5py("test_thread"))
    testSuite.addTest(test_h5py("test_process"))
    return testSuite


if __name__ == "__main__":
    import sys

    mysuite = main_test_suite()
    runner = unittest.TextTestRunner()
    if not runner.run(mysuite).wasSuccessful():
        sys.exit(1)

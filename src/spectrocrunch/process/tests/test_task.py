import unittest
from testfixtures import TempDirectory
from .. import utils


class test_task(unittest.TestCase):
    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()

    def _run_task(self, parameters, proc1, outputnxprocess=True):
        if proc1:
            previoustask = utils.nxpathtotask(proc1)
            self.assertTrue(previoustask.done)
        else:
            previoustask = None

        # Check run
        task2 = utils.create_task(dependencies=previoustask, **parameters)
        self.assertFalse(task2.done)
        task2.run()
        self.assertTrue(task2.done)
        if outputnxprocess:
            checksum = task2.checksum
            self.assertEqual(checksum, task2.output.checksum)
            self.assertTrue(task2.output.valid_checksum())
        proc2 = task2.output

        # Check re-run (same task instance)
        task2.run()
        proc3 = task2.output
        self.assertEqual(proc2, proc3)

        # Check re-run (new task instance)
        task3 = utils.create_task(dependencies=previoustask, **parameters)
        self.assertTrue(task3.done)
        task3.run()
        proc3 = task3.output
        self.assertEqual(proc2, proc3)

        # Check re-run (reconstructed task instance)
        if outputnxprocess:
            task4 = utils.nxpathtotask(proc2)
            self.assertEqual(type(task4), type(task2))
            self.assertTrue(task4.done)
            self.assertEqual(checksum, task4.checksum)
            self.assertEqual(checksum, task4.output.checksum)
            self.assertTrue(task4.output.valid_checksum())

            # Check re-run (new task instance)
            task4.run()
            proc3 = task4.output
            self.assertEqual(proc2, proc3)

        return proc2

    def _check_reproc(self, proc1, proc2):
        self.assertNotEqual(proc1, proc2)
        self.assertEqual(proc1.name.split(".")[-1], "1")
        self.assertEqual(proc2.name.split(".")[-1], "2")

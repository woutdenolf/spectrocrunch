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

        # Check run and re-run
        newtask = utils.create_task(dependencies=previoustask, **parameters)
        self.assertFalse(newtask.done)
        newtask.run()
        self.assertTrue(newtask.done)
        proc2 = newtask.output
        newtask.run()
        proc3 = newtask.output
        self.assertEqual(proc2, proc3)
        task = utils.create_task(dependencies=previoustask, **parameters)
        self.assertTrue(task.done)
        task.run()
        proc3 = task.output
        self.assertEqual(proc2, proc3)

        # Check reconstructed task from output
        if outputnxprocess:
            task = utils.nxpathtotask(proc2)
            self.assertEqual(type(task), type(newtask))
            self.assertTrue(task.done)
            task.run()
            proc3 = task.output
            self.assertEqual(proc2, proc3)

        return proc2

    def _check_reproc(self, proc1, proc2):
        self.assertNotEqual(proc1, proc2)
        self.assertEqual(proc1.name.split('.')[-1], '1')
        self.assertEqual(proc2.name.split('.')[-1], '2')

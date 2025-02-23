import logging
import multiprocessing
from ..utils import timing

logger = logging.getLogger(__name__)


def run_sequential(tasks, name=None):
    """
    Args:
        tasks(list(NXtask))
    Returns:
        bool
    """
    with timing.timeit_logger(logger, name=name):
        for task in tasks:
            task.run()
            if not task.done:
                return False
        return True


def run_task(task, result_queue):
    """
    Args:
        task(NXtask)
        result_queue(Queue)
    """
    try:
        task.run()
    finally:
        result_queue.put(task.done)


def run_parallel(tasks, name, nproc=2):
    """
    Args:
        tasks(list(NXtask))
    Returns:
        bool
    """
    with timing.timeit_logger(logger, name=name):
        result_queue = multiprocessing.Queue()
        results = []
        with multiprocessing.Pool(nproc) as pool:
            try:
                while tasks:
                    tasks_copy = list(tasks)
                    tasks = []
                    for task in tasks_copy:
                        if task.ready_to_run:
                            results.append(
                                pool.apply_async(run_task, task, result_queue)
                            )
                        else:
                            tasks.append(task)
                    if not result_queue.get():
                        for result in results:
                            result.wait()
                        return
            finally:
                for result in results:
                    result.wait()

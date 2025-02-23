from .basetask import TaskException
from ..utils import instance
from ..io import nxfs


def create_task(**parameters):
    method = parameters.get("method", None)
    if method == "crop":
        from .nxcrop import Task
    elif method == "replace":
        from .nxreplace import Task
    elif method == "minlog":
        from .nxminlog import Task
    elif method == "align":
        from .nxalign import Task
    elif method == "expression":
        from .nxexpression import Task
    elif method == "resample":
        from .nxresample import Task
    elif method == "pymca":
        from .nxpymca import Task
    elif method == "fullfield":
        from .nxfullfield import Task
    elif method == "xiaedftonx":
        from .nxxiaedf import Task
    elif method == "scenevis":
        from .scenevis import Task
    elif method == "xrfgeometry":
        from .nxqxrf import Task
    elif method == "blissmerge":
        from .nxblissmerge import Task
    elif method == "blissmcapre":
        from .nxblissmcapre import Task
    elif method == "blisspymca":
        from .nxblisspymca import Task
    elif method == "stack":
        from .nxstack import Task
    elif method == "tile":
        from .nxtile import Task
    else:
        Task = parameters.pop("_task", None)
        if Task is None:
            raise TaskException(
                "Unknown task defined by parameters {}".format(parameters)
            )
    return Task(**parameters)


def nxpathtotask(path):
    if instance.isstring(path):
        path = nxfs.factory(path)
    if path.is_nxclass("NXprocess"):
        parameters = path.config.read()
        if "method" not in parameters and "name" not in parameters:
            parameters["name"] = path.name
        from .nxprocesswrap import Task
    else:
        from .nxwrap import Task

        parameters = {"path": path}
    outputparent = path.parent
    dependencies = [path for path in path.dependencies]
    return create_task(
        dependencies=dependencies, outputparent=outputparent, _task=Task, **parameters
    )

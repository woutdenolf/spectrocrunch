# -*- coding: utf-8 -*-

import logging
from contextlib import contextmanager
import h5py
from ..io.utils import temporary_filename
from ..io.h5fs import prepare_h5data


logger = logging.getLogger(__name__)


@contextmanager
def open_uris(uris):
    sources = []
    files = []
    try:
        for path in uris:
            filename, path_in_file = path.split("::")
            f = h5py.File(filename, mode="r")
            files.append(f)
            sources.append(f[path_in_file])
        yield sources
    finally:
        for f in files:
            try:
                f.close()
            except Exception as e:
                logger.error("Closing '{}' : {}".format(f.filename, e))


def merge_h5datasets(
    dest, name, sources, shape_map, nscandim=1, no_stack_dimension=False
):
    """Merge datasets in a virtual dataset.

    :param h5py.Group dest:
    :param str name:
    :param list(h5py.Dataset) sources:
    :param dict shape_map:
    :param int nscandim: start index of the data dimensions
    """
    dset_shape = sources[0].shape
    dtype = sources[0].dtype
    fillvalue = sources[0].fillvalue

    lst = [dset.name for dset in sources if dset_shape != dset.shape]
    if lst:
        raise RuntimeError(f"Datasets do not have shape {dset_shape}: {lst}")

    if dset_shape:
        # VDS with reshaped scan dimensions
        data_shape = dset_shape[nscandim:]
        scan_shape = shape_map.get(dset_shape[:nscandim], dset_shape[:nscandim])
        # VDS does a C-order reshape while we need F-order
        scan_shape = scan_shape[::-1]
        if no_stack_dimension and len(sources) == 1:
            shape = scan_shape + data_shape
            layout = h5py.VirtualLayout(shape=shape, dtype=dtype)
            dset = sources[0]
            vsource = h5py.VirtualSource(
                dset.file.filename,
                dset.name,
                shape=dset.shape,
                dtype=dset.dtype,
            )
            layout[()] = vsource
        else:
            shape = (len(sources),) + scan_shape + data_shape
            layout = h5py.VirtualLayout(shape=shape, dtype=dtype)
            for i, dset in enumerate(sources):
                vsource = h5py.VirtualSource(
                    dset.file.filename,
                    dset.name,
                    shape=dset.shape,
                    dtype=dset.dtype,
                )
                layout[i] = vsource
        dest.create_virtual_dataset(name, layout, fillvalue=fillvalue)
    else:
        # Cannot make VDS of scalar datasets
        arr = [dset[()] for dset in sources]
        if isinstance(arr[0], str):
            if name == "end_time":
                arr = arr[-1]
            else:
                arr = arr[0]
        dest[name] = prepare_h5data(arr)


def merge_h5groups(
    dest_parent, dest_name, sources, shape_map, nscandim=1, no_stack_dimension=False
):
    """
    :param h5py.Group dest_parent:
    :param str dest_name:
    :param list(h5py.Group) sources:
    """
    source = sources[0]  # Assume HDF5-tree the same for all sources
    if dest_name in dest_parent:
        dest = dest_parent[dest_name]
    else:
        dest = dest_parent.create_group(dest_name)
        dest.attrs.update(source.attrs)

    for k in source:
        try:
            if k in dest:
                # Destination already exists: skip
                continue
            lnk = source.get(k, default=None, getlink=True)
            if isinstance(lnk, h5py.SoftLink):
                # Soft link: preserve the link path (destination changes)
                lnk_path = lnk.path
                if lnk_path.startswith("/"):
                    # TODO:
                    scan_name = "/".join(dest.name.split("/")[:2])
                    lnk_path = scan_name + "/" + "/".join(lnk_path.split("/")[2:])
                dest[k] = h5py.SoftLink(lnk_path)
                continue
            if isinstance(source[k], h5py.Group):
                # Group
                merge_h5groups(
                    dest,
                    k,
                    [s[k] for s in sources],
                    shape_map,
                    nscandim=nscandim,
                    no_stack_dimension=no_stack_dimension,
                )
            else:
                # Dataset
                merge_h5datasets(
                    dest,
                    k,
                    [s[k] for s in sources],
                    shape_map,
                    nscandim=nscandim,
                    no_stack_dimension=no_stack_dimension,
                )
        except Exception as e:
            logger.error(
                "Group '{}::{}' item '{}': {}".format(
                    dest.file.filename, dest.name, k, e
                )
            )
            raise


class MergedBlissScan(object):
    def __init__(self, scan_uris, filename=None, name=None):
        if name:
            self.name = name
        else:
            self.name = "merged"
        self.file = self._memfile_merged(scan_uris, filename=filename)

    @property
    def merged_scan(self):
        return self.file[self.name]

    def _memfile_merged(self, scan_uris, filename=None):
        if filename:
            kw = {}
        else:
            filename = temporary_filename(None, suffix=".h5")
            kw = {"driver": "core", "backing_store": False}
        h5file = h5py.File(filename, mode="a", **kw)
        try:
            with open_uris(scan_uris) as sources:
                merge_h5groups(
                    h5file, self.name, sources, self.shape_map(sources), nscandim=1
                )
            h5file.flush()
        except Exception:
            h5file.close()
            raise
        return h5file

    def shape_map(self, sources):
        return {}

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


class MergedBlissMesh(MergedBlissScan):
    def shape_map(self, sources):
        """Fast axis first"""
        if "title" not in sources[0]:
            return {}
        cmd = sources[0]["title"][()]
        for s in sources:
            assert cmd == s["title"][()]
        parts = cmd.split(" ")
        shape = int(parts[4]), int(parts[8]) + 1
        return {shape[0] * shape[1]: shape}

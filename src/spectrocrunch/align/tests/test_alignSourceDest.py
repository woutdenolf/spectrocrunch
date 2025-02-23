import unittest
from ..alignSource import alignSource
from ..alignDest import alignDest
from ..types import transformationType

from testfixtures import TempDirectory
from . import helper_teststack
import numpy as np
import os
import h5py
import fabio


class test_alignSource(unittest.TestCase):
    def setUp(self):
        self.dir = TempDirectory()

    def tearDown(self):
        self.dir.cleanup()

    def test_sourcedest(self):
        # Get test data
        listofstacks, offsets, stackdim = helper_teststack.data(
            transformationType.translation
        )
        nstacks = len(listofstacks)

        form = "stack%%0%dd" % int(np.floor(np.log10(nstacks)) + 1)
        names = [form % i for i in range(nstacks)]
        shape = listofstacks[0].shape
        nimages = np.take(shape, stackdim)
        imgsize = tuple(np.delete(shape, stackdim))
        dtype = listofstacks[0].dtype

        listofstacks_readonly = [s.copy() for s in listofstacks]
        listofstacks2 = [
            np.arange(8, dtype=dtype).reshape(2, 2, 2) for s in listofstacks
        ]
        listofstacks3 = [np.empty(s.shape, dtype=s.dtype) for s in listofstacks]
        listofstacks4 = [np.zeros(s.shape, dtype=s.dtype) for s in listofstacks]

        # Open HDF5 file
        h5name = "rawdata.h5"
        h5name2 = "rawdata2.h5"
        h5name3 = "rawdata3.h5"
        h5name4 = "rawdata4.h5"
        absh5name = os.path.join(self.dir.path, h5name)
        absh5name2 = os.path.join(self.dir.path, h5name2)
        absh5name3 = os.path.join(self.dir.path, h5name3)
        absh5name4 = os.path.join(self.dir.path, h5name4)
        with h5py.File(absh5name, mode="a") as h5f, h5py.File(
            absh5name2, mode="a"
        ) as h5f2, h5py.File(absh5name3, mode="a") as h5f3, h5py.File(
            absh5name4, mode="a"
        ) as h5f4:
            # Write h5 files
            for i in range(nstacks):
                h5f.create_dataset(names[i], data=listofstacks[i])
                h5f2.create_dataset(names[i], data=listofstacks2[i], maxshape=shape)
                h5f3.create_dataset(names[i], data=listofstacks3[i])
                h5f4.create_dataset(names[i], data=listofstacks4[i])
            self.dir.compare((h5name, h5name2, h5name3, h5name4))
            h5datasets = [h5f[name] for name in names]
            h5datasets2 = [h5f2[name] for name in names]
            h5datasets3 = [h5f3[name] for name in names]
            h5datasets4 = [h5f4[name] for name in names]

            # Write edf files
            form = "%%s_%%0%dd.edf" % int(np.floor(np.log10(nimages)) + 1)
            edfnames = [[form % (name, i) for i in range(nimages)] for name in names]
            absedfnames = [
                [
                    os.path.join(self.dir.path, "edfdir", form % (name, i))
                    for i in range(nimages)
                ]
                for name in names
            ]
            absedfnames2 = [
                [
                    os.path.join(self.dir.path, "edfdir2", form % (name, i))
                    for i in range(nimages)
                ]
                for name in names
            ]
            absedfnames3 = [
                [
                    os.path.join(self.dir.path, "edfdir3", form % (name, i))
                    for i in range(nimages)
                ]
                for name in names
            ]
            absedfnames_1 = [
                [
                    os.path.join(self.dir.path, "edfdir_1", form % (name, i))
                    for i in range(nimages)
                ]
                for name in names
            ]
            edfdir = os.path.join(self.dir.path, "edfdir")
            edfdir2 = os.path.join(self.dir.path, "edfdir2")
            edfdir3 = os.path.join(self.dir.path, "edfdir3")
            if not os.path.exists(edfdir):
                os.makedirs(edfdir)
            if not os.path.exists(edfdir2):
                os.makedirs(edfdir2)
            if not os.path.exists(edfdir3):
                os.makedirs(edfdir3)

            for i in range(nstacks):
                for j in range(nimages):
                    f = fabio.edfimage.edfimage(
                        data=np.take(listofstacks_readonly[i], j, axis=stackdim)
                    )
                    f.write(absedfnames[i][j])
                    f.write(absedfnames2[i][j])
            self.dir.compare(
                tuple([item for sublist in edfnames for item in sublist]), path="edfdir"
            )
            self.dir.compare(
                tuple([item for sublist in edfnames for item in sublist]),
                path="edfdir2",
            )

            # Check source construction
            source = [
                alignSource(absh5name, names, stackdim=stackdim),
                alignSource(
                    os.path.join(self.dir.path, "edfdir"), edfnames, stackdim=stackdim
                ),
                alignSource(listofstacks, None, stackdim=stackdim),
                alignSource(h5datasets, None, stackdim=stackdim),
            ]
            self.assertTrue(all(nimages == s.nimages for s in source))
            self.assertTrue(all(imgsize == s.imgsize for s in source))

            # Check source reading
            for i in range(nstacks):
                for j in range(nimages):
                    img = np.take(listofstacks_readonly[i], j, axis=stackdim)
                    for s in source:
                        np.testing.assert_array_equal(s.readimg(i, j), img)
                        np.testing.assert_array_equal(
                            s.readimgas(i, j, np.float32), img.astype(np.float32)
                        )

            # Test destination construction
            dest_exc = [
                alignDest(
                    listofstacks2, None, None, stackdim=stackdim
                ),  # edit datasets with wrong size
                alignDest(h5datasets2, None, None, stackdim=stackdim),
            ]  # edit datasets with wrong size

            dest = [
                alignDest(
                    absh5name, names, "", stackdim=stackdim
                ),  # adds an entry1 subdirectory
                alignDest(
                    absh5name2, names, "", stackdim=stackdim, overwrite=True
                ),  # overwrites the datasets
                # adds an entry1 subdirectory
                alignDest(absh5name3, names, "", stackdim=stackdim),
                alignDest(
                    absh5name3, names, ".aligned", stackdim=stackdim
                ),  # adds .aligned to names
                # adds .aligned to names in subdirectory entry1
                alignDest(absh5name3, names, ".aligned", stackdim=stackdim),
                # adds .aligned to names in subdirectory entry2
                alignDest(absh5name3, names, ".aligned", stackdim=stackdim),
                alignDest(
                    os.path.join(self.dir.path, "edfdir"),
                    names,
                    ".edf",
                    stackdim=stackdim,
                ),  # adds a ..._1 subdirectory
                alignDest(
                    os.path.join(self.dir.path, "edfdir2"),
                    names,
                    ".edf",
                    stackdim=stackdim,
                    overwrite=True,
                ),  # overwrites the files
                # creates files in the existing (empty) directory
                alignDest(
                    os.path.join(self.dir.path, "edfdir3"),
                    names,
                    ".edf",
                    stackdim=stackdim,
                ),
                alignDest(
                    listofstacks2, None, None, stackdim=stackdim, overwrite=True
                ),  # reform and edit np array
                alignDest(
                    listofstacks3, None, None, stackdim=stackdim
                ),  # edit np array
                alignDest(
                    h5datasets2, None, None, stackdim=stackdim, overwrite=True
                ),  # reform and edit datasets
                alignDest(h5datasets3, None, None, stackdim=stackdim),
            ]  # edit datasets

            # Check destination preparation
            for d in dest_exc:
                try:
                    d.prepare(nimages, imgsize, dtype)  # must throw an error
                    self.assertTrue(False)
                except ValueError:
                    pass
            for d in dest:
                d.prepare(nimages, imgsize, dtype)

            self.assertEqual(listofstacks2[0].shape, shape)

            # Check writing
            for i in range(nstacks):
                for j in range(nimages):
                    for d in dest:
                        img = np.take(listofstacks_readonly[i], j, axis=stackdim)
                        d.writeimg(img, i, j)

            # Check written files
            self.dir.compare(
                tuple([item for sublist in edfnames for item in sublist]), path="edfdir"
            )
            self.dir.compare(
                tuple([item for sublist in edfnames for item in sublist]),
                path="edfdir_1",
            )
            self.dir.compare(
                tuple([item for sublist in edfnames for item in sublist]),
                path="edfdir2",
            )
            self.dir.compare(
                tuple([item for sublist in edfnames for item in sublist]),
                path="edfdir3",
            )

            # Check written content
            h5datasets3_a = [h5f3["/entry1/" + name] for name in names]
            h5datasets3_b = [h5f3[name + ".aligned"] for name in names]
            h5datasets3_c = [h5f3["/entry1/" + name + ".aligned"] for name in names]
            h5datasets3_d = [h5f3["/entry2/" + name + ".aligned"] for name in names]
            for i in range(nstacks):
                for j in range(nimages):
                    img = np.take(listofstacks_readonly[i], j, axis=stackdim)
                    np.testing.assert_array_equal(
                        fabio.open(absedfnames[i][j]).data, img
                    )
                    np.testing.assert_array_equal(
                        fabio.open(absedfnames2[i][j]).data, img
                    )
                    np.testing.assert_array_equal(
                        fabio.open(absedfnames3[i][j]).data, img
                    )
                    np.testing.assert_array_equal(
                        fabio.open(absedfnames_1[i][j]).data, img
                    )
                np.testing.assert_array_equal(
                    listofstacks_readonly[i], listofstacks2[i]
                )
                np.testing.assert_array_equal(
                    listofstacks_readonly[i], listofstacks3[i]
                )
                np.testing.assert_array_equal(listofstacks_readonly[i], h5datasets[i])
                np.testing.assert_array_equal(listofstacks_readonly[i], h5datasets2[i])
                np.testing.assert_array_equal(listofstacks_readonly[i], h5datasets3[i])
                np.testing.assert_array_equal(
                    listofstacks_readonly[i], h5datasets3_a[i]
                )
                np.testing.assert_array_equal(
                    listofstacks_readonly[i], h5datasets3_b[i]
                )
                np.testing.assert_array_equal(
                    listofstacks_readonly[i], h5datasets3_c[i]
                )
                np.testing.assert_array_equal(
                    listofstacks_readonly[i], h5datasets3_d[i]
                )

                np.testing.assert_array_equal(
                    listofstacks_readonly[i] * 0, listofstacks4[i]
                )
                np.testing.assert_array_equal(
                    listofstacks_readonly[i] * 0, h5datasets4[i]
                )

    def test_teststack(self):
        listofstacks, offsets, stackdim = helper_teststack.data(
            transformationType.translation
        )
        self.assertIsInstance(listofstacks, list)
        self.assertIsInstance(listofstacks[0], np.ndarray)
        self.assertEqual(len(listofstacks[0].shape), 3)
        self.assertTrue(all(s.shape == listofstacks[0].shape for s in listofstacks))

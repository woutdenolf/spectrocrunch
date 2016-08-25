This project aims at processing spectrocopic imaging data (fluorescence and transmission).

Version
=======

    python setup.py version

Test
====

    python setup.py test

Run a single test:

    python -m spectrocrunch.align.tests.test_teststack

Dependencies on bare python 2.7
===============================
    
    pip install setuptools_scm
    pip install numpy
    pip install scipy
    pip install pylab 
    yum install hdf5-devel
    pip install h5py 
    pip install fabio 
    pip install testfixtures

    yum install cmake
    git clone https://github.com/kaspermarstal/SimpleElastix
    mkdir build
    cd build
    cmake ../SimpleElastix/SuperBuild
    make
    cd ./SimpleITK-build/Wrapping/PythonPackage
    python setup.py install

    pip install mako 
    pip install cffi
    pip install pyopencl
    wget https://github.com/woutdenolf/sift_pyocl/archive/master.zip
    unzip master
    cd sift_pyocl-master
    python setup.py install
        needed to change:
            add cl_buffer_region typedef in wrap_cl.h

    #pip install http://sourceforge.net/projects/pymca/files/pymca/PyMca5.1.1/


Installation
============

    python setup.py install

Use without installation
========================

To import modules from a package without installing the package, add the 
directory of the package to the PYTHONPATH environment variable or add this
to the top of your script
    import sys
    sys.path.insert(1,'/data/id21/inhouse/wout/dev/SpectroCrunch')

import as follows
    from spectrocrunch.align.alignElastix import alignElastix as align

Subpackages
===========

align
-----

    Aligning multiple image stacks with different alignment algorithms. One stack is the reference, the other stacks are aligned accordingly.

common
------

    Subpackage used by the other subpackages.

fullfield
---------

    Fullfield XAS data processing.

h5stacks
--------

    Data processing organized in a software independent hdf5 pipeline.

io
--

    Data I/O.

materials
---------

    Definition of compounds and mixtures with calculation of physical properties (database/calculation/simulation).

math
----

    Another subpackage used by the other subpackages, more specifically grouping all math.

process
-------

    This subpackage connects beamline specific code to the other subpackages.

visualization
-------------

    Plotting things.

xrf
---

    X-ray fluorescence data processing.



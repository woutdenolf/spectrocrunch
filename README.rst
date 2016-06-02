This project contains the following (sub)packages:

    ./
    |-- spectrocrunch/
    |   |-- align
    |   |   |-- tests
    |   |-- common
    |   |   |-- tests
    |   |-- io
    |   |   |-- tests
    |   |-- materials
    |   |   |-- tests
    |   |-- math
    |   |   |-- tests
    |   |-- process
    |   |   |-- tests
    |   |-- visualization
    |   |   |-- tests
    |   |-- xrf
    |   |   |-- tests
    |   |-- tests

Version
=======

    python setup.py version

Test
====

    python setup.py test

Run a single test:

    python -m spectrocrunch.align.tests.test_teststack

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

import as follows:
    from spectrocrunch.align.alignElastix import alignElastix as align

align
=====

    Aligning multiple image stacks with different alignment algorithms. One stack is the reference, the other stacks are aligned accordingly.

common
======

    Subpackage used by the other subpackages.

io
==

    Data I/O.

materials
=========

    Definition of compounds and mixtures with calculation of physical properties (database/calculation/simulation).

math
====

    Another subpackage used by the other subpackages, more specifically grouping all math.

process
=======

    This subpackage connects beamline specific code to the other subpackages.

visualization
=============

    Plotting things

xrf
===

    X-ray fluorescence data processing



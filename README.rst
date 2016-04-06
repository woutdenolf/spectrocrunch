This project contains the following (sub)packages:

    ./
    |-- spectrocrunch/
    |   |-- alignelastix
    |   |   |-- tests
    |   |-- tests


Use without installation
========================

To import modules from a package without installing the package, add the 
directory of the package to the PYTHONPATH environment variable and
import modules as follows:
    from spectrocrunch.align import align

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

alignelastix
============

    Align image stacks using SimpleElastix which is a Python binding for Elastix.

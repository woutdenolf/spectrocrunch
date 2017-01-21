SpectroCrunch: spectrocopic imaging library (XRF/XAS)
=====================================================

Getting started
---------------
Install:

.. code-block:: bash

    pip install spectrocrunch [--user]

Test:

.. code-block:: bash

    python -m spectrocrunch.test.test_all

Developers
----------
Main development website: https://github.com/woutdenolf/spectrocrunch

Distribution website: https://pypi.python.org/pypi/SpectroCrunch

Linux:   |Travis Status|
Windows: |Appveyor Status|

Guidelines for contributors and project managers can be found :ref:`dev-guide`

Use without installation
------------------------

.. code-block:: bash

    git clone https://github.com/woutdenolf/spectrocrunch
    cd spectrocrunch

To import modules from a package without installing the package, add the 
directory of the package to the PYTHONPATH environment variable or add this
to the top of your script

.. code-block::

    import sys
    sys.path.insert(1,'/data/id21/inhouse/wout/dev/SpectroCrunch')


import as follows:

.. code-block:: 

    from spectrocrunch.align.alignElastix import alignElastix as align

Subpackages
-----------

align
+++++

    Aligning multiple image stacks with different alignment algorithms. One stack is the reference, the other stacks are aligned accordingly.

common
++++++

    Subpackage used by the other subpackages.

fullfield
+++++++++

    Fullfield XAS data processing.

h5stacks
++++++++

    Data processing organized in a software independent hdf5 pipeline.

io
++

    Data I/O.

materials
+++++++++

    Definition of compounds and mixtures with calculation of physical properties (database/calculation/simulation).

math
++++

    Another subpackage used by the other subpackages, more specifically grouping all math.

process
+++++++

    This subpackage connects beamline specific code to the other subpackages.

visualization
+++++++++++++

    Plotting things.

xrf
+++

    X-ray fluorescence data processing.

.. |Travis Status| image:: https://travis-ci.org/woutdenolf/spectrocrunch.svg?branch=master
   :target: https://travis-ci.org/woutdenolf/spectrocrunch
.. |Appveyor Status| image:: https://ci.appveyor.com/api/projects/status/github/woutdenolf/spectrocrunch?svg=true
   :target: https://ci.appveyor.com/project/woutdenolf/spectrocrunch

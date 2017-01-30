SpectroCrunch: spectrocopic imaging library (XRF/XAS)
=====================================================

Getting started
---------------

Install dependencies:

.. code-block:: bash

    git clone https://github.com/woutdenolf/spectrocrunch
    . spectrocrunch/tools/prepare_installation.sh [-v 3]

Install from PyPi:

.. code-block:: bash

    pip install spectrocrunch [--user]

Install from source:

.. code-block:: bash

    python setup.py install [--user]

Test:

.. code-block:: bash

    python -m spectrocrunch.test.test_all

Developers
----------
Main development website: https://github.com/woutdenolf/spectrocrunch

Distribution website: https://pypi.python.org/pypi/SpectroCrunch

Guidelines for contributors and project managers can be found in tools/README.rst

+------------+-------------------------+-------------------------+
|            | Master                  | Develop                 |
+============+=========================+=========================+
| Linux      | |Travis Status Master|  | |Travis Status Develop| |
+------------+-------------------------+-------------------------+
| Windows    | |Appveyor Status|       |                         |
+------------+-------------------------+-------------------------+


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


.. |Travis Status Master| image:: https://travis-ci.org/woutdenolf/spectrocrunch.svg?branch=master
   :target: https://travis-ci.org/woutdenolf/spectrocrunch
.. |Travis Status Develop| image:: https://travis-ci.org/woutdenolf/spectrocrunch.svg?branch=develop
   :target: https://travis-ci.org/woutdenolf/spectrocrunch
.. |Appveyor Status| image:: https://ci.appveyor.com/api/projects/status/github/woutdenolf/spectrocrunch?svg=true
   :target: https://ci.appveyor.com/project/woutdenolf/spectrocrunch

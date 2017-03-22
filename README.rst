SpectroCrunch: spectrocopic imaging library (XRF/XAS)
=====================================================

Getting started
---------------

Install dependencies:

.. code-block:: bash

    git clone https://github.com/woutdenolf/spectrocrunch
    . spectrocrunch/tools/prepare_install-linux.sh

Install from PyPi:

.. code-block:: bash

    pip install spectrocrunch [--user]

Install from source:

.. code-block:: bash

    git clone https://github.com/woutdenolf/spectrocrunch Spectrocrunch
    cd Spectrocrunch
    pip install [--user] .

Test:

.. code-block:: bash

    python -m spectrocrunch.test.test_all

Developers
----------
Main development website: https://github.com/woutdenolf/spectrocrunch

Distribution website: https://pypi.python.org/pypi/SpectroCrunch

Pull requests must be done on the develop branch. More guidelines for contributors and project managers can be found in tools/README.rst

+------------+--------------------------+---------------------------+
|            | Master                   | Develop                   |
+============+==========================+===========================+
| Linux      | |Travis Status Master|   | |Travis Status Develop|   |
+------------+--------------------------+---------------------------+
| Windows    | |Appveyor Status Master| | |Appveyor Status Develop| |
+------------+--------------------------+---------------------------+


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
.. |Appveyor Status Master| image:: https://ci.appveyor.com/api/projects/status/github/woutdenolf/spectrocrunch?svg=true&branch=master
   :target: https://ci.appveyor.com/project/woutdenolf/spectrocrunch/branch/master
.. |Appveyor Status Develop| image:: https://ci.appveyor.com/api/projects/status/github/woutdenolf/spectrocrunch?svg=true&branch=develop
   :target: https://ci.appveyor.com/project/woutdenolf/spectrocrunch/branch/develop

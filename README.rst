SpectroCrunch: spectroscopic imaging library (XRF/XAS)
======================================================

Getting started
---------------

Install dependencies:

.. code-block:: bash

    git clone https://github.com/woutdenolf/spectrocrunch
    . spectrocrunch/tools/linux-install-deps.sh [-u]

Install from PyPi:

.. code-block:: bash

    pip install [--user] spectrocrunch

Install from source:

.. code-block:: bash

    git clone https://github.com/woutdenolf/spectrocrunch
    cd spectrocrunch
    pip install [--user] .

Test:

.. code-block:: bash

    # Test source:
    python setup.py test
    
    # Test installation:
    python -m spectrocrunch.tests.test_all
    
    # Test with options:
    python -m unittest -v spectrocrunch.tests.test_all.test_suite
    pytest --pyargs spectrocrunch
    nose2 -v spectrocrunch

Documentation:

 Release: http://pythonhosted.org/spectrocrunch

 Latest: http://spectrocrunch.readthedocs.io/en/latest/


Developers
----------

|Travis Status Master| |Appveyor Status Master|

Main development website: https://github.com/woutdenolf/spectrocrunch

Distribution website: https://pypi.python.org/pypi/spectrocrunch

Guidelines for contributors and project managers can be found in the `developers guide <https://github.com/woutdenolf/wdncrunch/blob/master/tools/README.rst/>`_.


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
    sys.path.insert(1,'.../spectrocrunch')


Import as follows:

.. code-block:: 

    from spectrocrunch.align.alignElastix import alignElastix as align


.. |Travis Status Master| image:: https://travis-ci.org/woutdenolf/spectrocrunch.svg?branch=master
   :target: https://travis-ci.org/woutdenolf/spectrocrunch
.. |Appveyor Status Master| image:: https://ci.appveyor.com/api/projects/status/github/woutdenolf/spectrocrunch?svg=true&branch=master
   :target: https://ci.appveyor.com/project/woutdenolf/spectrocrunch/branch/master

Guide for developers
====================

Versioning
----------

Versioning is done on the master branch::

  git checkout master
  git merge develop --no-commit 
  #... change the version according to the rules below ...
  git add .
  git commit -m "Bump version to 1.2.3-beta4"

Semantic versioning is followed (http://semver.org/)::

  MAJOR.MINOR.MICRO.SERIAL

  SERIAL: bump when changes not to the code
  MICRO : bump when bug fix is done
              when bumping SERIAL == 15
  MINOR : bump when API changes backwards compatible
              when new functionality is added
              when bumping MICRO == 15
  MAJOR : bump when API changes not backwards compatible
 
  Always reset the lower numbers to 0, except for SERIAL which starts at 1.

  dev   : not tested
  alpha : begin testing
  beta  : feature complete
  rc    : test complete
  final : stable version


Releasing
---------

A release can be created on github with the tag MAJOR.MINOR.MICRO

TODO: upon releasing on github, trigger ci with pypi deployment

From the source
---------------

.. code-block:: bash

    git clone https://github.com/woutdenolf/spectrocrunch

    . spectrocrunch/tools/prepare_installation.sh [-v 3]
    if [[ $? == 0 ]]; then echo "OK"; else echo "NOT OK"; fi

    cd spectrocrunch
    python setup.py version
    python setup.py test
    python -m spectrocrunch.align.tests.test_teststack

    python setup.py build
    python setup.py install [--user]
    # OR
    pip install . [--user]
    
Manual Deployment
-----------------

Add PyPi credentials file ~/.pypirc (chmod 600):

.. code-block:: bash

    [distutils]
    index-servers =
      pypi
      pypitest

    [pypi]
    repository=https://pypi.python.org/pypi
    username=...
    password=...

    [pypitest]
    repository=https://testpypi.python.org/pypi
    username=...
    password=...


Register project:

.. code-block:: bash

    python setup.py register -r pypi
    python setup.py register -r pypitest

Deploy:

.. code-block:: bash

    # on linux
    python setup.py sdist bdist_wheel upload -r pypi
    # on windows
    python setup.py bdist_msi upload -r pypi
    
Help
----

.. code-block:: bash

    python setup.py --help-commands
    python setup.py sdist --help-formats
    python setup.py bdist --help-formats

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

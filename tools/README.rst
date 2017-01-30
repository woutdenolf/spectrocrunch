Guide for developers
====================

Git branching model
-------------------

`This <http://nvie.com/posts/a-successful-git-branching-model/>`_ branching model is followed:

* Github branches: master (RELEV=final), develop (RELEV=dev)

* Local branches: master (RELEV=final), develop (RELEV=dev), feat-\*, fix-\*, hotfix-\*

.. code-block:: bash

  git clone https://github.com/woutdenolf/spectrocrunch Spectrocrunch
  git config --global push.followTags true

* Feature/fix branch:

.. code-block:: bash

  # Start working on the feature
  git checkout -b feat-something develop

  # Commit changes ...

  # Merge feature into develop
  git checkout develop
  git merge --no-ff feat-something
  git branch -d feat-something

  # Publish
  git push origin develop

* Release branch:

.. code-block:: bash

  # Start releasing
  git checkout -b release-1.2 develop

  # change RELEV from "dev" to "alpha"
  # Possible change in SERIAL, MICRO (bug fixes) and RELEV (testing progress)

  # Version to be released
  echo `python -c "from _version import version;print(\"v{}\".format(version));"`

  # CHANGELOG.rst: add release

  # Bump the version
  git add .
  git commit -m "Bump version to 1.2.3"

  # Merge release in master and develop
  git checkout master
  git merge --no-ff release-1.2
  git tag -s v1.2.3 -m "Version 1.2.3"

  git checkout develop
  git merge --no-ff release-1.2

  git branch -d release-1.2

  # Publish
  git push origin develop
  git push origin master

* Hotfix branch:

.. code-block:: bash

  # Master tag is 1.2.3
  git checkout -b hotfix-1.2.4 master

  # Possible change in SERIAL, MICRO (bug fixes) and RELEV (testing progress)
  # Finally RELEV=final

  # CHANGELOG.rst: add release

  git add .

  # Check current version
  echo `python -c "from _version import version;print(\"v{}\".format(version));"`

  # Bump the version
  git commit -m "Bump version to 1.2.4"

  # Merge release in master and develop
  git checkout master
  git merge --no-ff hotfix-1.2.4
  git tag -s v1.2.4 -m "Version 1.2.4"
  git checkout develop
  git merge --no-ff hotfix-1.2.4

  git branch -d hotfix-1.2.4

Versioning
----------

`Semantic versioning <http://semver.org/>`_ is followed::

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

Create a release on github based on this tag

  Title: Release of version MAJOR.MINOR.MICRO

  Body: Copy from CHANGELOG

   
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


From the source
---------------

.. code-block:: bash

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

Guide for developers
====================

Versioning
----------

Versioning is done on the master branch::

  git checkout master
  git merge develop --no-commit 
  #... change the version according to the rules below ...
  git add .
  git commit -m "bump version to 1.2.3-beta4"

Semantic versioning is followed (http://semver.org/)::

  MAJOR.MINOR.MICRO.SERIAL

  SERIAL: bump when changes not to the code
  MICRO : bump when bug fix is done
              when SERIAL == 15
  MINOR : bump when API changes backwards compatible
              when new functionality is added
              when MICRO == 15
  MAJOR : bump when API changes not backwards compatible
 
  Always reset the lower numbers to 0.

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

    git clone https://github.com/woutdenolf/spectrocrunch
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

    # define pypi and pypitest in .pypirc (chmod 600 .pypirc)
    python setup.py register -r pypi
    # OR
    python setup.py register -r testpypi

    # on linux
    python setup.py sdist  bdist_wheel upload -r pypi
    # on windows
    python setup.py bdist_msi upload -r pypi
    
Help
----

    python setup.py --help-commands
    python setup.py sdist --help-formats
    python setup.py bdist --help-formats

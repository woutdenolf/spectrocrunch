Guide for developers
====================

Configure git
-------------

.. code-block:: bash

    git config --global user.name githubname
    git config --global user.email user@domain
    git config --global user.signingkey YOURMASTERKEYID

(See `signing  <localrefsigning_>`_)


Contribute
----------

Assuming you forked spectrocrunch on github, then your fork will be referred to as "origin" and the repository you forked from will be referred to as "upstream".

* Clone your fork locally:

.. code-block:: bash

  git clone https://github.com/forkuser/spectrocrunch
  git remote add upstream https://github.com/woutdenolf/spectrocrunch

* Add a feature:

.. code-block:: bash

  # Branch of the upstream master
  git fetch upstream
  git checkout upstream/master
  git branch feat-something
  git checkout feat-something

  # Stay up to date with upstream
  git branch --set-upstream-to=upstream/master
  git pull

  # Commit changes ...
  git commit -m "..."

  # Push to origin
  git push origin feat-something

* Create a new pull request with

  base fork: woutdenolf/spectrocrunch (upstream)

  base: master

  head fork: forkuser/spectrocrunch (origin)

  compare: feat-something

* Keep your master up to date:

.. code-block:: bash
  
  git checkout master
  git pull upstream master (== git fetch upstream; git merge upstream/master)
  git push origin master

* Clean up your repository:

.. code-block:: bash
  
  git fetch -p upstream


Increase the version
--------------------

1. Get the master

.. code-block:: bash
  
  git checkout master
  git pull upstream master

2. Update version in _version.py and update CHANGELOG.rst (see `versioning  <localrefversion_>`_)

.. code-block:: bash
  
  echo `python -c "from _version import version;print(\"v{}\".format(version));"`

3. Check whether the version is in releasable state (see `check releasable state  <localrefreleasable_>`_)

4. Commit and tag new version

.. code-block:: bash
  
  git add .
  git commit -m "Bump version to 1.2.3"
  git tag -s v1.2.3 -m "Version 1.2.3"
  push origin
  push origin v1.2.3

5. Create a new pull request with

  base fork: woutdenolf/spectrocrunch (upstream)

  base: master

  head fork: forkuser/spectrocrunch (origin)

  compare: v1.2.3


.. _localrefreleasable:

Check releasable state
----------------------

1. Create a clean `sandbox <localrefsandbox_>`_ and make a fresh git clone

2. Release directory

.. code-block:: bash
  
  export RELEASEDIR=...
  export VERSION=`python -c "from _version import strictversion as version;print(\"{}\".format(version));"`
  rm -r ${RELEASEDIR}
  mkdir -p ${RELEASEDIR}/dist

3. Build the source tarball

.. code-block:: bash
  
  python setup.py clean sdist
  cp dist/spectrocrunch-${VERSION}.tar.gz ${RELEASEDIR}/dist

4. Test the source

.. code-block:: bash
  
  tar zxvf ${RELEASEDIR}/dist/spectrocrunch-${VERSION}.tar.gz
  cd spectrocrunch-${VERSION}
  pip install .
  python -m spectrocrunch.tests.test_all
  
5. Release the docs

.. code-block:: bash
  
  python setup.py clean build_doc
  pip uninstall -y spectrocrunch
  cd build/sphinx/html
  zip -r ${RELEASEDIR}/html_doc.zip .
  cd ../../..

6. Inspect the docs

.. code-block:: bash
  
  firefox build/sphinx/html/index.html

7. Build the wheels (do this on different platforms)

.. code-block:: bash
  
  python setup.py clean bdist_wheel
  cp dist/spectrocrunch-${VERSION}-py2.py3-none-any.whl ${RELEASEDIR}/dist

8. Test the wheels

.. code-block:: bash
  
  pip install ${RELEASEDIR}/dist/spectrocrunch-${VERSION}-py2.py3-none-any.whl
  python -m spectrocrunch.tests.test_all
  pip uninstall -y spectrocrunch

9. Delete the `sandbox  <localrefsandbox_>`_


Release a version
-----------------

1. Get the version to be released

.. code-block:: bash
  
  git checkout master
  git pull upstream master
  git checkout v1.2.3

2. Check whether the version is in releasable state (see `check releasable state  <localrefreleasable_>`_). This should have been done when creating the release but best to double check.

3. Create a release on github based on the tag

  Title: Release of version MAJOR.MINOR.MICRO

  Body: Copy from CHANGELOG

4. Deploy code (see `pypi setup  <localrefdeployment_>`_)

.. code-block:: bash

  twine upload -r pypitest --sign ${RELEASEDIR}/*
  twine upload -r pypi --sign ${RELEASEDIR}/*

5. Deploy documentation

.. code-block:: bash

  https://testpypi.python.org/pypi?%3Aaction=pkg_edit&name=spectrocrunch
  http://pypi.python.org/pypi?%3Aaction=pkg_edit&name=spectrocrunch


.. _localrefdeployment:

Deployment
----------

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

Register project (already done):

.. code-block:: bash

  twine register -r pypi dist/*.whl
  twine register -r pypitest dist/*.whl


.. _localrefsandbox:

Sandbox
-------

* Using virtualenv

.. code-block:: bash

  virtualenv test1.2.3
  cd test1.2.3
  source bin/activate

* Using pyenv

Installation and activation

.. code-block:: bash

  export PYTHON_CONFIGURE_OPTS="--enable-shared"
  export PYENV_ROOT="~/.pyenv"
  if [[ ! -d $PYENV_ROOT ]]; then
    git clone https://github.com/pyenv/pyenv.git ${PYENV_ROOT}
    git clone https://github.com/pyenv/pyenv-virtualenv.git ${PYENV_ROOT}/plugins/pyenv-virtualenv
  fi
  export PATH="$PYENV_ROOT/bin:$PATH"
  eval "$(pyenv init -)"
  eval "$(pyenv virtualenv-init -)"

Manage python versions

  pyenv install 2.7.13
  pyenv uninstall 2.7.13

  pyenv local 2.7.13 (in this directory)
  pyenv shell 2.7.13 (in this shell)
  pyenv shell --unset

  pyenv version
  pyenv versions

Manage virtualenvs

  pyenv virtualenv 2.7.13 myenvname
  pyenv activate myenvname
  pyenv deactivate
  pyenv uninstall myenvname
  pyenv virtualenvs

.. _localrefversion:

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
 
  Always reset the lower numbers to 0.

  dev   : not tested
  alpha : begin testing
  beta  : feature complete
  rc    : test complete
  final : stable version


Install external dependencies
-----------------------------

.. code-block:: bash

    . spectrocrunch/tools/prepare_install-linux.sh [-v 3]
    if [[ $? == 0 ]]; then echo "OK"; else echo "NOT OK"; fi


Help
----

.. code-block:: bash

    python setup.py --help-commands
    python setup.py sdist --help-formats
    python setup.py bdist --help-formats


.. _localrefsigning:

Signing
-------

Generate PGP keypair:

.. code-block:: bash

    while true; do ls -R / &>/dev/null; sleep 1; done &
    gpg --gen-key

Generate a revocation certificate:

.. code-block:: bash

    gpg --output revoke.asc --gen-revoke YOURMASTERKEYID
    shred --remove revoke.asc

Publish public key:

.. code-block:: bash

    gpg --keyserver pgp.mit.edu --send-keys YOURMASTERKEYID

Share public key:

.. code-block:: bash

    gpg --armor --export YOURMASTERKEYID
    (or look it up in pgp.mit.edu)

Revoke PGP key:

.. code-block:: bash

    gpg --keyserver pgp.mit.edu --recv-keys YOURMASTERKEYID
    gpg --import revoke.asc
    gpg --keyserver pgp.mit.edu --send-keys YOURMASTERKEYID

Share private PGP key:

.. code-block:: bash

    gpg --export-secret-key -a | ssh user@host gpg --import -

Show all keys:

.. code-block:: bash

    gpg --list-keys


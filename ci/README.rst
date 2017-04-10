Guide for continuous integration
================================

1. `Travis CI <localreftravis_>`_ (linux)

2. `AppVeyor <localrefappveyor_>`_ (windows)


.. _localreftravis:

Travis CI
---------

Pre-build dependencies
++++++++++++++++++++++

1. Install Ubuntu 12.04.5 LTS: http://releases.ubuntu.com/12.04/

.. code-block:: bash

  adduser travis
  usermod -aG sudo travis

2. Install Python 2.7.12 from source:

.. code-block:: bash

  sudo -s
  apt-get update
  apt-get install libbz2-dev libsqlite3-dev libreadline-dev zlib1g-dev libncurses5-dev libssl-dev libgdbm-dev libssl-dev openssl tk-dev

  wget https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz
  tar -xvzf Python-2.7.12.tgz
  rm Python-2.7.12.tgz
  cd Python-2.7.12

  ./configure --prefix=/opt/python/2.7.12 --enable-shared LDFLAGS=-Wl,-rpath=/opt/python/2.7.12/lib
  make -j2
  make install

  cd ..
  rm -rf Python-2.7.12

  export PATH="/opt/python/2.7.12/bin/:$PATH"

3. Install pip:

.. code-block:: bash

  python -m ensurepip
  pip install --upgrade pip
  pip install virtualenv

4. Create virtual env:

.. code-block:: bash

  cd /home/travis
  virtualenv python2.7.12
  source python2.7.12/bin/activate

5. Install dependencies:

.. code-block:: bash

  apt-get install git
  git clone https://github.com/woutdenolf/spectrocrunch.git
  . spectrocrunch/tools/prepare_install-linux.sh -u

Accept when:

.. code-block:: bash

  Python version: 2.7.12 
  Python location: /home/travis/virtualenv/python2.7.12/bin/python 
  Python include: /opt/python/2.7.12/include/python2.7 
  Python library: /opt/python/2.7.12/lib/libpython2.7.so 
  Pip: 9.0.1 from /home/travis/virtualenv/python2.7.12/lib/python2.7/site-packages (python 2.7) 

6. Create pre-build and upload:

.. code-block:: bash

  tar -czf spectrocrunch.travis.python2.7.tgz 2.7/simpleelastix 2.7/cmake
  curl --upload-file spectrocrunch.travis.python2.7.tgz https://transfer.sh/spectrocrunch.travis.python2.7.tgz

.. _localrefappveyor:

AppVeyor
--------


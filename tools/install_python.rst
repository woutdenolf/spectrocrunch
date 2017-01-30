Use pyenv
---------

.. code-block:: bash

    pyenv install 3.4.0
    pyenv versions # show the versions you can choose for, including the system
    pyenv version # show current python version

    # python for this user
    pyenv global version
    pyenv global
    
    # python in directory and all subdirectories
    pyenv local version
    pyenv local --unset
    pyenv local

    # python in current shell
    pyenv shell version 
    pyenv shell --unset
    pyenv shell

Install system-wide Python locally (NOT recommended, use pyenv instead)
-----------------------------------------------------------------------

.. code-block:: bash

    apt-get install libssl-dev openssl tk-dev
    cd /opt
    wget https://www.python.org/ftp/python/3.4.6/Python-3.4.6.tgz
    tar -xvzf Python-3.4.6.tgz
    cd Python-3.4.6

    ./configure --enable-optimizations --enable-shared
    make -j2
    make install

    # python will complain about missing library
    # solution 1:
    #   LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
    # solution 2:
    #   configure again without --enable-shared

    ./configure --enable-optimizations
    make -j2
    make install


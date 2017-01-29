
Install Python 3 on Debian 7 (tested with 3.4.6 and 3.5.3):

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




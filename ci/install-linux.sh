#!/bin/bash

cwd=$(pwd)
echo $cwd

echo "Installing curl ..."
$version = curl_version();
$ssl_supported= ($version['features'] & CURL_VERSION_SSL);
echo $version
echo $version['features']
echo $ssl_supported

#./configure --with-ssl
#make
#sudo make install
#cd $cwd

echo "Installing cmake ..."
wget --no-check-certificate http://www.cmake.org/files/v3.2/cmake-3.2.2.tar.gz
tar xf cmake-3.2.2.tar.gz
cd cmake-3.2.2
./configure
make
sudo make install
cd $cwd

echo "Installing SimpleElastix ..."
mkdir simpleelastix
cd simpleelastix
git clone -b master https://github.com/kaspermarstal/SimpleElastix
mkdir build
cd build
cmake ../SimpleElastix/SuperBuild
make
cd ./SimpleITK-build/Wrapping/Python/Packaging
python setup.py install
cd $cwd


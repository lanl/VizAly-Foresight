#!/bin/bash

# SZ
echo "Building SZ ... "

git clone https://github.com/disheng222/SZ.git
cd SZ
git checkout v2.1.4.2
mkdir install
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..

echo "Building SZ done!"

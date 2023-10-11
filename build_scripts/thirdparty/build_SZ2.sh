#!/bin/bash

echo "Building SZ 2 ... "

git clone https://github.com/disheng222/SZ.git
cd SZ
git checkout v2.1.12.5
mkdir install
mkdir build
cd build

cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..

echo "Building SZ 2.1.12.5 done!"

#!/bin/bash

echo "Building SZ 3 ... "

git clone https://github.com/szcompressor/SZ3.git
cd SZ3
git checkout v3.1.7
mkdir install
mkdir build
cd build

cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..

echo "Building SZ 3.1.7 done!"

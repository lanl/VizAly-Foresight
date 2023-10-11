#!/bin/bash

echo "Building SPERR ... "


git clone https://github.com/NCAR/SPERR.git
cd SPERR/

git checkout v0.7.1
mkdir build
mkdir install
cd build

cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..

echo "Building SPERR SZ 0.7.1 done!"
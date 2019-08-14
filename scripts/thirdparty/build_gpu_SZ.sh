#!/bin/bash

# SZ
echo "Building GPU version of SZ ... "

wget https://lanl.github.io/VizAly-Foresight/compressors/SZ-GPU-generic.tar.gz
tar -zxvf SZ-GPU-generic.tar.gz
cd SZ-generic
mkdir install
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..

echo "Building GPU SZ done!"

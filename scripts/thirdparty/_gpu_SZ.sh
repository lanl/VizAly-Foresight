#!/bin/bash

# SZ
echo "Building GPU version of SZ ... "

# git clone https://github.com/robertu94/SZ.git SZ-generic
# git checkout generic
wget https://lanl.github.io/VizAly-Foresight/compressors/SZ-GPU-generic.tar.gz
mkdir SZ-generic-GPU
tar -zxvf SZ-GPU-generic.tar.gz --directory SZ-generic-GPU
cd SZ-generic-GPU
mkdir install
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..

echo "Building GPU SZ done!"

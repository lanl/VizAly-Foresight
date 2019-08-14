#!/bin/bash

# SZ
echo "Building GPU version of SZ ... "

wget https://github.com/lanl/VizAly-Foresight/blob/gh-pages/compressors/SZ-GPU-generic.zip
unzip ./SZ-GPU-generic.zip
rm -f ./SZ-GPU-generic.zip
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
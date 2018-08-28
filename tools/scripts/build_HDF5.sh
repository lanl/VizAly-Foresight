#!/bin/bash

# HDF5
echo "Building HDF5 ... "

git https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git
cd hdf5
git checkout hdf5-1_10_3
mkdir install
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..

echo "Building HDF5 done!"
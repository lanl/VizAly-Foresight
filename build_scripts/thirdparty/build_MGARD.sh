#!/bin/bash

# MGARD
echo "Building MGARD ... "

git clone https://github.com/CODARcode/MGARD.git
cd MGARD

git checkout 0f6cdf9b59e837547e3298be7187de8df376bb18

mkdir install
mkdir build

cd build

echo "Building MGARD"
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j

echo "Installing MGARD"
make install
cd ..
cd ..

echo "Building MGARD done!"
#!/bin/bash

# MGARD
echo "Building MGARD ... "

git clone https://github.com/CODARcode/MGARD.git
cd MGARD

mkdir install
mkdir build

cd build

echo "Installing MGARD"

cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..

echo "Building MGARD done!"

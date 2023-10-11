#!/bin/bash

echo "Building MGARD ... "

git clone https://github.com/CODARcode/MGARD.git
cd MGARD
git checkout 1.5.2

mkdir install
mkdir build

cd build


cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j

echo "Installing MGARD"
make install
cd ..
cd ..

echo "Building MGARD 1.5.2 done!"
#!/bin/bash

# SZ
echo "Building AMReX ... "

git clone https://github.com/AMReX-Codes/amrex
cd amrex

mkdir install
mkdir build
cd build

cmake .. -DCMAKE_INSTALL_PREFIX=../install -DBUILD_SHARED_LIBS=ON -DENABLE_AMRDATA=ON -DENABLE_OMP=ON
make -j
make install

cd ..

cd ..

echo "Building AMReX done!"

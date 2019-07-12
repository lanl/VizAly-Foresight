#!/bin/bash

# BigCrunch
echo "Building fpzip ... "

wget https://computation.llnl.gov/projects/floating-point-compression/download/fpzip-1.2.0.tar.gz
tar -zxvf fpzip-1.2.0.tar.gz
cd fpzip-1.2.0
cd src
make CXXFLAGS+='-fPIC' -j
cd ..
cd ..

echo "Building fpzip done!"
#!/bin/bash

# BigCrunch
echo "Building fpzip ... "

#wget https://computing.llnl.gov/projects/floating-point-compression/download/fpzip-1.2.0.tar.gz
wget https://lanl.github.io/VizAly-Foresight/compressors/fpzip-1.2.0.tar.gz

tar -zxvf fpzip-1.2.0.tar.gz
cd fpzip-1.2.0
cd src
make CXXFLAGS+='-fPIC' -j
cd ..
cd ..

echo "Building fpzip done!"
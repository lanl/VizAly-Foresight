#!/bin/bash

# BigCrunch
echo "Building BigCrunch ... "

git clone https://github.com/lanl/VizAly-BigCrunch.git
cd VizAly-BigCrunch/
mkdir install
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..

echo "Building BigCrunch done!"
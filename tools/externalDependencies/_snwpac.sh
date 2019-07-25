#!/bin/bash

# SnowPack
echo "Building snwpac ... "

git clone https://github.com/jpulidojr/VizAly-SNWPAC.git
cd VizAly-SNWPAC/
#git checkout v0.3 #Use master for now
mkdir install
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install -DSNWPAC_BUILD_STATIC:BOOL=false
make -j 
make install
cd ..
cd ..

echo "Building snwpac done!"

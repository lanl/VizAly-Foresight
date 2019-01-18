#!/bin/bash

# LossyWave
echo "Building LossyWave ... "

git clone https://github.com/lanl/VizAly-LossyWave.git
cd VizAly-LossyWave/
git checkout v0.2
mkdir install
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install -DLOSSYWAVE_BUILD_STATIC:BOOL=false
make -j 
make install
cd ..
cd ..

echo "Building LossyWave done!"

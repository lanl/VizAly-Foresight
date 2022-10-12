#!/bin/bash

# BigCrunch
echo "Building zfp ... "

git clone https://github.com/LLNL/zfp.git
cd zfp/
git checkout 1.0.0
mkdir install
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..

echo "Building ZFP done!"

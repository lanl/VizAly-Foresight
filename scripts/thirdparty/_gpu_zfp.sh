#!/bin/bash

# BigCrunch
echo "Building cuda_zfp ... "

git clone https://github.com/LLNL/zfp.git gpu_zfp
cd gpu_zfp/
git checkout 0.5.4
mkdir install
mkdir build
cd build
cmake -DZFP_WITH_CUDA=ON .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..

echo "Building CUDA_ZFP done!"
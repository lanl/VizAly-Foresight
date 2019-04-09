#!/bin/bash 

projectPath=$(pwd)

# Get parameter 
if [ "$#" -ne 1 ]; then
    buildDir=build
else
  buildDir=$1
fi

# Create build directory 
mkdir $buildDir
cd $buildDir

# build 
cmake ../src -DCMAKE_C_FLAGS=-dynamic -DCMAKE_CXX_FLAGS=-dynamic \
    -DCBENCH_ENABLE_BLOSC=OFF \
    -DBLOSC_INCLUDE_PATH=$projectPath/ExternalDependencies/c-blosc/install/include \
    -DBLOSC_LIBRARY=$projectPath/ExternalDependencies/c-blosc/install/lib/libblosc.so \
    -DCBENCH_ENABLE_SZ=ON \
    -DSZ_INCLUDE_PATH=$projectPath/ExternalDependencies/SZ/sz/include \
    -DSZ_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libSZ.a \
    -DZLIB_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libzlib.a \
    -DZSTD_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libzstd.a \
    -DCBENCH_ENABLE_LOSSY_WAVE=ON \
    -DLOSSYWAVE_INCLUDE_PATH=$projectPath/ExternalDependencies/VizAly-LossyWave/install/include \
    -DLOSSYWAVE_LIBRARY=$projectPath/ExternalDependencies/VizAly-LossyWave/install/lib/liblossywave.a \
    -DLOSSYWAVE_LZ4_LIBRARY=$projectPath/ExternalDependencies/VizAly-LossyWave/build/3rdparty/lz4/lz4-external/src/lz4-external/lib/liblz4.a \
    -DCBENCH_ENABLE_ZFP=ON \
    -DZFP_INCLUDE_PATH=$projectPath/ExternalDependencies/zfp/install/include \
    -DZFP_LIBRARY=$projectPath/ExternalDependencies/zfp/install/lib64/libzfp.a \
    -DCBENCH_ENABLE_FPZIP=ON \
    -DFPZIP_INCLUDE_PATH=$projectPath/ExternalDependencies/fpzip-1.2.0/inc/ \
    -DFPZIP_LIBRARY=$projectPath/ExternalDependencies/fpzip-1.2.0/lib/libfpzip.a \
    -DCBENCH_ENABLE_ISABELA=ON \
    -DISABELA_INCLUDE_PATH=$projectPath/ExternalDependencies/ISABELA-compress-0.2.1/include \
    -DISABELA_LIBRARY=$projectPath/ExternalDependencies/ISABELA-compress-0.2.1/lib/libisabela.a \
    -DCBENCH_ENABLE_BIG_CRUNCH=OFF
make -j

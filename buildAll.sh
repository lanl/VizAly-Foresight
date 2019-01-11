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
cmake ../src \
    -DBLOSC_INCLUDE_PATH=$projectPath/ExternalDependencies/c-blosc/install/include \
	-DBLOSC_LIBRARY=$projectPath/ExternalDependencies/c-blosc/install/lib/libblosc.so \
	-DCBENCH_ENABLE_SZ=ON \
	-DSZ_INCLUDE_PATH=$projectPath/ExternalDependencies/SZ/sz/include \
	-DSZ_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libSZ.so \
	-DZLIB_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libzlib.so \
	-DZSTD_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libzstd.so \
	-DCBENCH_ENABLE_BIG_CRUNCH=ON \
	-DBIGCRUNCH_INCLUDE_PATH=$projectPath/ExternalDependencies/VizAly-BigCrunch/install/include \
	-DBIGCRUNCH_LIBRARY=$projectPath/ExternalDependencies/VizAly-BigCrunch/install/lib/libbigcrunch.so \
	-DCBENCH_ENABLE_LOSSY_WAVE=ON \
	-DLOSSYWAVE_INCLUDE_PATH=$projectPath/ExternalDependencies/VizAly-LossyWave/install/include \
	-DLOSSYWAVE_LIBRARY=$projectPath/ExternalDependencies/VizAly-LossyWave/install/lib/liblossywave.so \
	-DCBENCH_ENABLE_ZFP=ON \
	-DZFP_INCLUDE_PATH=$projectPath/ExternalDependencies/zfp/install/include \
	-DZFP_LIBRARY=$projectPath/ExternalDependencies/zfp/install/lib/libzfp.so \
	-DCBENCH_ENABLE_NYX_LOADER=ON \
	-DHDF5_DIR=$projectPath/ExternalDependencies/hdf5/install/share/cmake/hdf5
make -j


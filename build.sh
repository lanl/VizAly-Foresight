#!/bin/bash 

cd ExternalDependencies
source buildExternal.sh
cd ..
 
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
cmake ../src -DBLOSC_INCLUDE_PATH=../ExternalDependencies/c-blosc/install/include -DBLOSC_LIBRARY=../ExternalDependencies/c-blosc/install/lib/libblosc.so -DSZ_INCLUDE_PATH=../ExternalDependencies/SZ/sz/include -DSZ_LIBRARY=../ExternalDependencies/SZ/install/lib/libSZ.so -DZLIB_LIBRARY=../ExternalDependencies/SZ/install/lib/libzlib.so -DZSTD_LIBRARY=../ExternalDependencies/SZ/install/lib/libzstd.so
make -j16 

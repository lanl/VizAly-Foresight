#!/bin/bash

# HDF5
echo "Building HDF5 ... "

echo $PLATFORM

if [ "$PLATFORM" = "travis" ]; then
	echo "Travis: Using docker HDF5 build"
else
	git clone https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git
	cd hdf5
	git checkout hdf5-1_10_3
	mkdir install
	mkdir build
	cd build
	#cmake .. -DHDF5_BUILD_CPP_LIB:BOOL=false -DHDF5_ENABLE_PARALLEL:BOOL=true -DCMAKE_C_FLAGS=-dynamic -DCMAKE_CXX_FLAGS=-dynamic -DCMAKE_INSTALL_PREFIX=../install - on CRAY
	cmake .. -DHDF5_BUILD_CPP_LIB:BOOL=false -DHDF5_ENABLE_PARALLEL:BOOL=true -DCMAKE_INSTALL_PREFIX=../install
	make -j
	make install
	cd ..
	cd ..
fi

echo "Building HDF5 done!"

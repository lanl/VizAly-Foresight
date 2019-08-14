#!/bin/bash

# Builds dependencies
currentDir=$(pwd)

for file in build_gpu_*.sh;
do
    pushd "$currentDir"  # to guard against build failure
    	source $file;
    popd
done

source _HDF5.sh;
source build_fpzip.sh;

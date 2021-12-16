#!/bin/bash

# Builds dependencies
currentDir=$(pwd)

for file in *_*.sh;
do
	if [[ ${file} != *"gpu_"* ]]  && [[ ${file} != *"__"* ]] ;then
    	pushd "$currentDir"  # to guard against build failure
    		source $file;
    	popd
	fi
done

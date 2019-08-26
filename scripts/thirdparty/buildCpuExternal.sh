#!/bin/bash

# Builds dependencies
currentDir=$(pwd)

for file in *_*.sh;
do
	if [[ ${testmystring} != *"_gpu_"* ]];then
    	pushd "$currentDir"  # to guard against build failure
    		source $file;
    	popd
	fi
done

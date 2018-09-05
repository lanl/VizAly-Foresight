#!/bin/bash

# Builds dependencies
currentDir=$(pwd)

for file in build_*.sh;
do
    pushd "$currentDir"  # to guard against build failure
    	source $file;
    popd
done

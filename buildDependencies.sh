#!/bin/bash 
projectPath=$(pwd)

# create folder for dependencies
mkdir ExternalDependencies
cd ExternalDependencies

# copy the scripts into that folder and run 
cp tools/scripts/*.sh .
source buildExternal.sh
cd ..
 

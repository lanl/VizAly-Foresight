#!/bin/bash 
 
# Get parameter 
if [ "$#" -ne 1 ]; then 
    buildDir=build 
else 
  buildDir=$1 
fi 

#create install directory
mkdir install

# Create build directory 
mkdir $buildDir 
cd $buildDir 
 
# build 
cmake ../src 
make -j16 
make install
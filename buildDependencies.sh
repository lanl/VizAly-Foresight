#!/bin/bash 
projectPath=$(pwd)


# Default arguments
buildType="Default"


# Process arguments
for arg in "$@"; do
	if [ $arg = "-all" ]; then
		buildType="all"

	fi

	if [ $arg = "-h" ]; then
		echo "Arguments: nothing for standard build, \"-all\" for all, and -h for help"
		return
	fi
done


# create folder for dependencies
mkdir ExternalDependencies
cd ExternalDependencies


# copy the scripts into that folder and run 
cp ../tools/externalDependencies/*.sh .
if [ $buildType = "all" ]; then
	echo "building all dependencies ..."
	source buildAllExternal.sh
else
	source buildExternal.sh
fi
cd ..

#
# Build Arguments
# -all : build all dependencies
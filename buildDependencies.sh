#!/bin/bash 
projectPath=$(pwd)


# Default arguments
buildType="Default"
externalDependenciesPath=$projectPath"/ExternalDependencies"


# Process arguments
index=0
for ((i=1; i<=$#; i++ )); do
	arg=${!i}

	if [ $arg = "-h" ]; then
		echo ""
		echo "Build Arguments:"
		echo "  --path <folder> : expernal dependencies path"
		echo ""
		echo " --all: build all dependencies"
		echo ""
		return
	fi

	if [ $arg = "-all" ]; then
		buildType="all"
	fi

	if [ $arg = "--path" ]; then
		index=i
		index=$((index+1))
		externalDependenciesPath=${!index}
	fi
done


# create folder for dependencies
mkdir $externalDependenciesPath
pushd $externalDependenciesPath


# copy the scripts into that folder and run 
cp $projectPath/tools/externalDependencies/*.sh .
if [ $buildType = "all" ]; then
	echo "building all dependencies ..."
	source buildAllExternal.sh
else
	source buildExternal.sh
fi
popd

#
# Build Arguments
# -all   : build all dependencies
# --path : build dependency path

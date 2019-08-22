#!/bin/bash 
projectPath=$(pwd)


# Default arguments
buildType="Default"
externalDependenciesPath=$projectPath"/ExternalDependencies"


# Process arguments
index=0
for ((i=1; i<=$#; i++ )); do
	arg=${!i}

	if [ $arg = "-h" ] || [ $arg = "--help" ]; then
		echo ""
		echo "Build Arguments:"
		echo " --all: build all possible dependencies"
		echo " --cpu: build dependencies for a CPU build"
		echo " --gpu: build dependencies for a GPU build"
		echo " --path <folder> : external dependencies path"
		echo ""
		return
	fi

	if [ $arg = "--all" ]; then
		buildType="all"
	fi

	if [ $arg = "--cpu" ]; then
		buildType="cpu"
	fi

	if [ $arg = "--gpu" ]; then
		buildType="gpu"
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
cp $projectPath/scripts/thirdparty/*.sh .
if [ $buildType = "all" ]; then
	echo "building all possible dependencies ..."
	source buildAllExternal.sh
elif [ $buildType = "cpu" ]; then
	echo "building all CPU dependencies ..."
	source buildCpuExternal.sh
elif [ $buildType = "gpu" ]; then
	echo "building all GPU dependencies ..."
	source buildGpuExternal.sh
else
	source buildExternal.sh
fi
popd

#
# Build Arguments:
# --all  : build all dependencies
# --cpu  : build CPU dependencies
# --gpu  : build GPU dependencies
# --path : build dependency path

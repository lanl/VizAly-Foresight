#!/bin/bash 
projectPath=$(pwd)


# Default arguments
buildType="cpu"
externalDependenciesPath=$projectPath"/ExternalDependencies"


# Process arguments
index=0
for ((i=1; i<=$#; i++ )); do
	arg=${!i}

	if [ $arg = "-h" ] || [ $arg = "--help" ]; then
		echo ""
		echo "Build Arguments:"
		echo " --cpu: build dependencies for a CPU build (default)"
		echo " --gpu: build dependencies for a GPU build"
		echo " --path <folder> : external dependencies path (default folder is ExternalDependencies)"
		echo ""
		return
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
cp $projectPath/build_scripts/thirdparty/*.sh .
if [ $buildType = "cpu" ]; then
	echo "building all CPU dependencies ..."
	source buildCpuExternal.sh
elif [ $buildType = "gpu" ]; then
	echo "building all GPU dependencies ..."
	source buildGpuExternal.sh
else
	echo "building default dependencies ..."
	source buildCpuExternal.sh
fi
popd

#
# Build Arguments:
# --cpu  : build CPU dependencies (default)
# --gpu  : build GPU dependencies
# --path : build dependency path (default folder is ExternalDependencies)

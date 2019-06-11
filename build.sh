#!/bin/bash 
projectPath=$(pwd)

# Default arguments
buildDir="build"
buildType="Debug"
buildOpt="Default"
buildPlatform="Default"
opts=""


# Process arguments
index=0
for ((i=1; i<=$#; i++ )); do
	arg=${!i}

	#help
	if [ $arg = "-h" ]; then

		echo ""
		echo "Build Arguments:"
		echo "  --path <folder> : build exec in that folder"
		echo ""
		echo "  -all : build with all options on"
		echo "  -min : build with minimal options on"
		echo "  -hacc: build for HACC only; no HDF5"
		echo ""
		echo "  -cori   : build for cori at nersc"
		echo "  -cooley : build for cooley"
		echo ""
		echo "  -release: use release mode instead of debug"
		echo ""
		echo "  --cleanAll <folder>: remove build folder and ExternalDependencies"
		echo "  --clean <folder>: remove build folder"
		echo ""
		echo "  No arguments: build path is build, and type is Debug"
		return
	fi

	# build location
	if [ $arg = "--path" ]; then
		index=i
		index=$((index+1))
		buildDir=${!index}
	fi

	if [ $arg = "--clean" ]; then
		if  [ -z "$2" ]; then
			echo "removing build"
			rm -rf build
		else
			index=i
			index=$((index+1))
			buildDir=${!index}

			echo "removing folder" $buildDir
			rm -rf $buildDir
		fi
		return
	fi

	# clean
	if [ $arg = "--cleanAll" ]; then
		next_arg=$((arg+1))
		if  [ -z "$2" ]; then
			echo "removing build"
			rm -rf build
		else
			index=i
			index=$((index+1))
			buildDir=${!index}

			echo "removing folder" $buildDir
			rm -rf $buildDir
		fi

		rm -rf ExternalDependencies
		return
	fi

	# platform
	if [ $arg = "-cori" ]; then
		buildPlatform="cori"
	fi

	# type
	if [ $arg = "-release" ]; then
		buildType="Release"
	fi

	


	# options
	if [ $arg = "-all" ]; then
		buildOpt="all"
	fi


	if [ $arg = "-hacc" ]; then
		buildOpt="cori"
	fi

	if [ $arg = "-min" ]; then
		buildOpt="minimal"
	fi

	if [ $arg = "-osx" ]; then
		buildOpt="osx"
	fi

	index=$((index+1))
done
 


# Platform
if [ $buildPlatform = "cori" ]; then
	export CPATH=/usr/common/software/hdf5-parallel/1.10.1/gnu/include:$CPATH
	export LD_LIBRARY_PATH=/usr/common/software/hdf5-parallel/1.10.1/gnu/lib:$LD_LIBRARY_PATH

	opt = "-DCMAKE_C_FLAGS=-dynamic -DCMAKE_CXX_FLAGS=-dynamic "

	mv CBench/CMakeLists.txt CBench/CMakeLists.txt.old
	mv scripts/Cori.CMakeLists.txt CBench/CMakeLists.txt
fi


# Create build directory 
mkdir $buildDir 
cd $buildDir 
 


# build
if [ $buildOpt = "min" ]; then
	echo "Default with minimal options ..."
	cmake ../CBench \
		-DCMAKE_BUILD_TYPE=$buildType

elif [ $buildOpt = "Default" ]; then
	echo "Default build on ..."
	cmake ../CBench $opt\
	    -DCBENCH_ENABLE_NYX_LOADER=ON \
		-DHDF5_DIR=$projectPath/ExternalDependencies/hdf5/install/share/cmake/hdf5 \
	    -DCBENCH_ENABLE_BLOSC=ON \
	    -DBLOSC_INCLUDE_PATH=$projectPath/ExternalDependencies/c-blosc/install/include \
		-DBLOSC_LIBRARY=$projectPath/ExternalDependencies/c-blosc/install/lib/libblosc.so \
		-DCBENCH_ENABLE_SZ=ON \
		-DSZ_INCLUDE_PATH=$projectPath/ExternalDependencies/SZ/sz/include \
		-DSZ_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libSZ.so \
		-DZLIB_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libzlib.so \
		-DZSTD_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libzstd.so \
		-DCBENCH_ENABLE_ZFP=ON \
		-DZFP_INCLUDE_PATH=$projectPath/ExternalDependencies/zfp/install/include \
		-DZFP_LIBRARY=$projectPath/ExternalDependencies/zfp/install/lib64/libzfp.so \
		-DCBENCH_ENABLE_FPZIP=ON \
		-DFPZIP_INCLUDE_PATH=$projectPath/ExternalDependencies/fpzip-1.2.0/inc/ \
		-DFPZIP_LIBRARY=$projectPath/ExternalDependencies/fpzip-1.2.0/lib/libfpzip.a \
		-DCBENCH_ENABLE_ISABELA=ON \
		-DISABELA_INCLUDE_PATH=$projectPath/ExternalDependencies/ISABELA-compress-0.2.1/include \
		-DISABELA_LIBRARY=$projectPath/ExternalDependencies/ISABELA-compress-0.2.1/lib/libisabela.a \
		-DCMAKE_BUILD_TYPE=$buildType

elif [ $buildOpt = "hacc" ]; then
	echo "Building for HACC, no hdf5 ..."

	cmake ../CBench $opt\
	    -DCBENCH_ENABLE_BLOSC=ON \
	    -DBLOSC_INCLUDE_PATH=$projectPath/ExternalDependencies/c-blosc/install/include \
		-DBLOSC_LIBRARY=$projectPath/ExternalDependencies/c-blosc/install/lib/libblosc.so \
		-DCBENCH_ENABLE_SZ=ON \
		-DSZ_INCLUDE_PATH=$projectPath/ExternalDependencies/SZ/sz/include \
		-DSZ_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libSZ.so \
		-DZLIB_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libzlib.so \
		-DZSTD_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libzstd.so \
		-DCBENCH_ENABLE_ZFP=ON \
		-DZFP_INCLUDE_PATH=$projectPath/ExternalDependencies/zfp/install/include \
		-DZFP_LIBRARY=$projectPath/ExternalDependencies/zfp/install/lib64/libzfp.so \
		-DCBENCH_ENABLE_FPZIP=ON \
		-DFPZIP_INCLUDE_PATH=$projectPath/ExternalDependencies/fpzip-1.2.0/inc/ \
		-DFPZIP_LIBRARY=$projectPath/ExternalDependencies/fpzip-1.2.0/lib/libfpzip.a \
		-DCBENCH_ENABLE_ISABELA=ON \
		-DISABELA_INCLUDE_PATH=$projectPath/ExternalDependencies/ISABELA-compress-0.2.1/include \
		-DISABELA_LIBRARY=$projectPath/ExternalDependencies/ISABELA-compress-0.2.1/lib/libisabela.a \
		-DCMAKE_BUILD_TYPE=$buildType

elif [ $buildOpt = "all" ]; then
	echo "Building with all dependencies ..."

	cmake ../CBench $opt\
	    -DCBENCH_ENABLE_NYX_LOADER=ON \
		-DHDF5_DIR=$projectPath/ExternalDependencies/hdf5/install/share/cmake/hdf5 \
	    -DCBENCH_ENABLE_BLOSC=ON \
	    -DBLOSC_INCLUDE_PATH=$projectPath/ExternalDependencies/c-blosc/install/include \
		-DBLOSC_LIBRARY=$projectPath/ExternalDependencies/c-blosc/install/lib/libblosc.so \
		-DCBENCH_ENABLE_SZ=ON \
		-DSZ_INCLUDE_PATH=$projectPath/ExternalDependencies/SZ/sz/include \
		-DSZ_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libSZ.so \
		-DZLIB_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libzlib.so \
		-DZSTD_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libzstd.so \
		-DCBENCH_ENABLE_ZFP=ON \
		-DZFP_INCLUDE_PATH=$projectPath/ExternalDependencies/zfp/install/include \
		-DZFP_LIBRARY=$projectPath/ExternalDependencies/zfp/install/lib64/libzfp.so \
		-DCBENCH_ENABLE_FPZIP=ON \
		-DFPZIP_INCLUDE_PATH=$projectPath/ExternalDependencies/fpzip-1.2.0/inc/ \
		-DFPZIP_LIBRARY=$projectPath/ExternalDependencies/fpzip-1.2.0/lib/libfpzip.a \
		-DCBENCH_ENABLE_ISABELA=ON \
		-DISABELA_INCLUDE_PATH=$projectPath/ExternalDependencies/ISABELA-compress-0.2.1/include \
		-DISABELA_LIBRARY=$projectPath/ExternalDependencies/ISABELA-compress-0.2.1/lib/libisabela.a \
		-DCBENCH_ENABLE_MGARD=ON \
		-DMGARD_INCLUDE_PATH=$projectPath/ExternalDependencies/MGARD/install/include \
		-DMGARD_LIBRARY=$projectPath/ExternalDependencies/MGARD/install/lib/libmgard.so \
		-DCBENCH_ENABLE_LOSSY_WAVE=ON \
		-DLOSSYWAVE_INCLUDE_PATH=$projectPath/ExternalDependencies/VizAly-LossyWave/install/include \
		-DLOSSYWAVE_LIBRARY=$projectPath/ExternalDependencies/VizAly-LossyWave/install/lib/liblossywave.so \
		-DLOSSYWAVE_LZ4_LIBRARY=$projectPath/ExternalDependencies/VizAly-LossyWave/build/3rdparty/lz4/lz4-external/src/lz4-external/lib/liblz4.a \
		-DCBENCH_ENABLE_BIG_CRUNCH=ON \
		-DBIGCRUNCH_INCLUDE_PATH=$projectPath/ExternalDependencies/VizAly-BigCrunch/install/include \
		-DBIGCRUNCH_LIBRARY=$projectPath/ExternalDependencies/VizAly-BigCrunch/install/lib/libbigcrunch.so \
		-DCMAKE_BUILD_TYPE=$buildType

elif [ $buildOpt = "osx" ]; then
	echo "Building for osx ..."

	# OSX for Hoby
	cmake ../CBench \
	    -DCMAKE_BUILD_TYPE=Debug \
	    -DCBENCH_ENABLE_BLOSC=ON \
	    -DBLOSC_INCLUDE_PATH=$projectPath/ExternalDependencies/c-blosc/install/include \
		-DBLOSC_LIBRARY=$projectPath/ExternalDependencies/c-blosc/install/lib/libblosc.a \
		-DCBENCH_ENABLE_SZ=ON \
		-DSZ_INCLUDE_PATH=$projectPath/ExternalDependencies/SZ/sz/include \
		-DSZ_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libSZ.dylib \
		-DZLIB_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libzlib.dylib \
		-DZSTD_LIBRARY=$projectPath/ExternalDependencies/SZ/install/lib/libzstd.dylib \
		-DCBENCH_ENABLE_BIG_CRUNCH=ON \
		-DBIGCRUNCH_INCLUDE_PATH=$projectPath/ExternalDependencies/VizAly-BigCrunch/install/include \
		-DBIGCRUNCH_LIBRARY=$projectPath/ExternalDependencies/VizAly-BigCrunch/install/lib/libbigcrunch.dylib \
		-DCBENCH_ENABLE_LOSSY_WAVE=ON \
		-DLOSSYWAVE_INCLUDE_PATH=$projectPath/ExternalDependencies/VizAly-LossyWave/install/include \
		-DLOSSYWAVE_LIBRARY=$projectPath/ExternalDependencies/VizAly-LossyWave/install/lib/liblossywave.dylib \
		-DLOSSYWAVE_LZ4_LIBRARY=$projectPath/ExternalDependencies/VizAly-LossyWave/build/3rdparty/lz4/lz4-external/src/lz4-external/lib/liblz4.a \
		-DCBENCH_ENABLE_MGARD=OFF \
		-DMGARD_INCLUDE_PATH=$projectPath/ExternalDependencies/MGARD/install/include \
		-DMGARD_LIBRARY=$projectPath/ExternalDependencies/MGARD/install/lib/libmgard.dylib \
		-DCBENCH_ENABLE_ZFP=ON \
		-DZFP_INCLUDE_PATH=$projectPath/ExternalDependencies/zfp/install/include \
		-DZFP_LIBRARY=$projectPath/ExternalDependencies/zfp/install/lib/libzfp.dylib \
		-DCBENCH_ENABLE_NYX_LOADER=ON \
		-DHDF5_DIR=$projectPath/ExternalDependencies/hdf5/install/share/cmake/hdf5 \
		-DCBENCH_ENABLE_FPZIP=ON \
		-DFPZIP_INCLUDE_PATH=$projectPath/ExternalDependencies/fpzip-1.2.0/inc/ \
		-DFPZIP_LIBRARY=$projectPath/ExternalDependencies/fpzip-1.2.0/lib/libfpzip.a \
		-DCBENCH_ENABLE_ISABELA=ON \
		-DISABELA_INCLUDE_PATH=$projectPath/ExternalDependencies/ISABELA-compress-0.2.1/include \
		-DISABELA_LIBRARY=$projectPath/ExternalDependencies/ISABELA-compress-0.2.1/lib/libisabela.a
else
	echo "Build type " $buildOpt " not supported yet!"
fi

make -j


#
# Build Arguments
# --path <foder name> : build exec in that folder

# -cooley : build for cooley
# -cori : build for cori at nersc

# -all : build with all options on
# -min : build with minimal options on
# -hacc : build for HACC only; no HDF5

# -release: use release mode instead of debug

# --cleanAll <folder>: remove build folder and ExternalDependencies
# --clean <folder>: remove build folder

# No argument is default:
# 	build path is build
# 	build type is Debug
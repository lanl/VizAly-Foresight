#!/bin/bash 
projectPath=$(pwd)

# Default arguments
buildDir="build"
buildType="Debug"
buildOpt="Default"
buildPlatform="Default"
opts=""
externalDependencies=$projectPath"/ExternalDependencies"


# Process arguments
index=0
for ((i=1; i<=$#; i++ )); do
	arg=${!i}

	#help
	if [ $arg = "-h" ] || [ $arg = "--help" ]; then

		echo ""
		echo "Build Arguments:"
		echo "  --path <folder> : build exec in that folder"
		echo ""
		echo " --externalDependencies <folder>: path where external dependencies was built"
		echo ""
		echo "  --all : build with all options on"
		echo "  --min : build with minimal options on"
		echo "  --hacc: build for HACC only; no HDF5"
		echo ""
		echo "  --cori   : build for cori at nersc"
		echo "  --cooley : build for cooley"
		echo ""
		echo "  --release: use release mode instead of debug"
		echo ""
		echo "  --clean <folder>: remove build folder"
		echo "  --cleanAll <folder>: remove build folder and ExternalDependencies"
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

	# external dependencies location
	if [ $arg = "--externalDependencies" ]; then
		index=i
		index=$((index+1))
		externalDependencies=${!index}
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
	if [ $arg = "--cori" ]; then
		buildPlatform="cori"
	fi

	# type
	if [ $arg = "--release" ]; then
		buildType="Release"
	fi

	


	# options
	if [ $arg = "--all" ]; then
		buildOpt="all"
	fi

	if [ $arg = "--travis" ]; then
		buildOpt="travis"
	fi

	if [ $arg = "--hacc" ]; then
		buildOpt="hacc"
	fi

	if [ $arg = "--min" ]; then
		buildOpt="minimal"
	fi

	if [ $arg = "--osx" ]; then
		buildOpt="osx"
	fi

	if [ $arg = "--gpu" ]; then
		buildOpt="gpu"
	fi

	index=$((index+1))
done
 


# Platform
if [ $buildPlatform = "cori" ]; then
	echo ""
	export CPATH=/opt/cray/pe/hdf5-parallel/1.10.2.0/GNU/8.2/include:$CPATH
	export LD_LIBRARY_PATH=/opt/cray/pe/hdf5-parallel/1.10.2.0/GNU/8.2/lib:$LD_LIBRARY_PATH

	#opt="-DCMAKE_C_FLAGS=-dynamic -DCMAKE_CXX_FLAGS=-dynamic -DCBENCH_PLATFORM=Cori "
	opt="-DCBENCH_PLATFORM=Cori"
else
    opt="-DCBENCH_PLATFORM=Default "
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
	if [ $buildPlatform = "cori" ]; then
		echo "..Building on Cori.."
		cmake ../CBench $opt\
			-DCBENCH_ENABLE_NYX_LOADER=ON \
			-DHDF5_DIR=$externalDependencies/hdf5/install/share/cmake/hdf5 \
			-DCBENCH_ENABLE_BLOSC=ON \
			-DBLOSC_INCLUDE_PATH=$externalDependencies/c-blosc/install/include \
			-DBLOSC_LIBRARY=$externalDependencies/c-blosc/install/lib/libblosc.so \
			-DCBENCH_ENABLE_SZ=ON \
			-DSZ_INCLUDE_PATH=$externalDependencies/SZ/sz/include \
			-DSZ_LIBRARY=$externalDependencies/SZ/install/lib/libSZ.a \
			-DZLIB_LIBRARY=$externalDependencies/SZ/install/lib/libzlib.a \
			-DZSTD_LIBRARY=$externalDependencies/SZ/install/lib/libzstd.a \
			-DCBENCH_ENABLE_ZFP=ON \
			-DZFP_INCLUDE_PATH=$externalDependencies/zfp/install/include \
			-DZFP_LIBRARY=$externalDependencies/zfp/install/lib64/libzfp.a \
			-DCBENCH_ENABLE_FPZIP=ON \
			-DFPZIP_INCLUDE_PATH=$externalDependencies/fpzip-1.2.0/inc/ \
			-DFPZIP_LIBRARY=$externalDependencies/fpzip-1.2.0/lib/libfpzip.a \
			-DCBENCH_ENABLE_ISABELA=ON \
			-DISABELA_INCLUDE_PATH=$externalDependencies/ISABELA-compress-0.2.1/include \
			-DISABELA_LIBRARY=$externalDependencies/ISABELA-compress-0.2.1/lib/libisabela.a \
			-DCMAKE_BUILD_TYPE=$buildType
	else
		cmake ../CBench $opt\
			-DCBENCH_ENABLE_NYX_LOADER=ON \
			-DHDF5_DIR=$externalDependencies/hdf5/install/share/cmake/hdf5 \
			-DCBENCH_ENABLE_BLOSC=ON \
			-DBLOSC_INCLUDE_PATH=$externalDependencies/c-blosc/install/include \
			-DBLOSC_LIBRARY=$externalDependencies/c-blosc/install/lib/libblosc.so \
			-DCBENCH_ENABLE_SZ=ON \
			-DSZ_INCLUDE_PATH=$externalDependencies/SZ/sz/include \
			-DSZ_LIBRARY=$externalDependencies/SZ/install/lib/libSZ.so \
			-DZLIB_LIBRARY=$externalDependencies/SZ/install/lib/libzlib.so \
			-DZSTD_LIBRARY=$externalDependencies/SZ/install/lib/libzstd.so \
			-DCBENCH_ENABLE_ZFP=ON \
			-DZFP_INCLUDE_PATH=$externalDependencies/zfp/install/include \
			-DZFP_LIBRARY=$externalDependencies/zfp/install/lib64/libzfp.so \
			-DCBENCH_ENABLE_FPZIP=ON \
			-DFPZIP_INCLUDE_PATH=$externalDependencies/fpzip-1.2.0/inc/ \
			-DFPZIP_LIBRARY=$externalDependencies/fpzip-1.2.0/lib/libfpzip.a \
			-DCBENCH_ENABLE_ISABELA=ON \
			-DISABELA_INCLUDE_PATH=$externalDependencies/ISABELA-compress-0.2.1/include \
			-DISABELA_LIBRARY=$externalDependencies/ISABELA-compress-0.2.1/lib/libisabela.a \
			-DCMAKE_BUILD_TYPE=$buildType
	fi

elif [ $buildOpt = "hacc" ]; then
	echo "Building for HACC, no hdf5 ..."

	cmake ../CBench $opt\
		-DCBENCH_ENABLE_BLOSC=ON \
		-DBLOSC_INCLUDE_PATH=$externalDependencies/c-blosc/install/include \
		-DBLOSC_LIBRARY=$externalDependencies/c-blosc/install/lib/libblosc.so \
		-DCBENCH_ENABLE_SZ=ON \
		-DSZ_INCLUDE_PATH=$externalDependencies/SZ/sz/include \
		-DSZ_LIBRARY=$externalDependencies/SZ/install/lib/libSZ.so \
		-DZLIB_LIBRARY=$externalDependencies/SZ/install/lib/libzlib.so \
		-DZSTD_LIBRARY=$externalDependencies/SZ/install/lib/libzstd.so \
		-DCBENCH_ENABLE_ZFP=ON \
		-DZFP_INCLUDE_PATH=$externalDependencies/zfp/install/include \
		-DZFP_LIBRARY=$externalDependencies/zfp/install/lib64/libzfp.so \
		-DCBENCH_ENABLE_FPZIP=ON \
		-DFPZIP_INCLUDE_PATH=$externalDependencies/fpzip-1.2.0/inc/ \
		-DFPZIP_LIBRARY=$externalDependencies/fpzip-1.2.0/lib/libfpzip.a \
		-DCBENCH_ENABLE_ISABELA=ON \
		-DISABELA_INCLUDE_PATH=$externalDependencies/ISABELA-compress-0.2.1/include \
		-DISABELA_LIBRARY=$externalDependencies/ISABELA-compress-0.2.1/lib/libisabela.a \
		-DCMAKE_BUILD_TYPE=$buildType

elif [ $buildOpt = "all" ]; then
	echo "Building with all dependencies ..."

	if [ "$PLATFORM" = "travis" ]; then
	  echo "Travis: Using internal HDF5 build"
	  cmake ../CBench $opt\
		-DCBENCH_ENABLE_NYX_LOADER=ON \
		-DCBENCH_ENABLE_BLOSC=ON \
		-DBLOSC_INCLUDE_PATH=$externalDependencies/c-blosc/install/include \
		-DBLOSC_LIBRARY=$externalDependencies/c-blosc/install/lib/libblosc.so \
		-DCBENCH_ENABLE_SZ=ON \
		-DSZ_INCLUDE_PATH=$externalDependencies/SZ/sz/include \
		-DSZ_LIBRARY=$externalDependencies/SZ/install/lib/libSZ.so \
		-DZLIB_LIBRARY=$externalDependencies/SZ/install/lib/libzlib.so \
		-DZSTD_LIBRARY=$externalDependencies/SZ/install/lib/libzstd.so \
		-DCBENCH_ENABLE_ZFP=ON \
		-DZFP_INCLUDE_PATH=$externalDependencies/zfp/install/include \
		-DZFP_LIBRARY=$externalDependencies/zfp/install/lib64/libzfp.so \
		-DCBENCH_ENABLE_FPZIP=ON \
		-DFPZIP_INCLUDE_PATH=$externalDependencies/fpzip-1.2.0/inc/ \
		-DFPZIP_LIBRARY=$externalDependencies/fpzip-1.2.0/lib/libfpzip.a \
		-DCBENCH_ENABLE_ISABELA=ON \
		-DISABELA_INCLUDE_PATH=$externalDependencies/ISABELA-compress-0.2.1/include \
		-DISABELA_LIBRARY=$externalDependencies/ISABELA-compress-0.2.1/lib/libisabela.a \
		-DCBENCH_ENABLE_MGARD=ON \
		-DMGARD_INCLUDE_PATH=$externalDependencies/MGARD/install/include \
		-DMGARD_LIBRARY=$externalDependencies/MGARD/install/lib/libmgard.so \
		-DCBENCH_ENABLE_LOSSY_WAVE=ON \
		-DLOSSYWAVE_INCLUDE_PATH=$externalDependencies/VizAly-LossyWave/install/include \
		-DLOSSYWAVE_LIBRARY=$externalDependencies/VizAly-LossyWave/install/lib/liblossywave.so \
		-DLOSSYWAVE_LZ4_LIBRARY=$externalDependencies/VizAly-LossyWave/build/3rdparty/lz4/lz4-external/src/lz4-external/lib/liblz4.a \
		-DCBENCH_ENABLE_BIG_CRUNCH=ON \
		-DBIGCRUNCH_INCLUDE_PATH=$externalDependencies/VizAly-BigCrunch/install/include \
		-DBIGCRUNCH_LIBRARY=$externalDependencies/VizAly-BigCrunch/install/lib/libbigcrunch.so \
		-DCMAKE_BUILD_TYPE=$buildType
	else
	  cmake ../CBench $opt\
		-DCBENCH_ENABLE_NYX_LOADER=ON \
		-DHDF5_DIR=$externalDependencies/hdf5/install/share/cmake/hdf5 \
		-DCBENCH_ENABLE_BLOSC=ON \
		-DBLOSC_INCLUDE_PATH=$externalDependencies/c-blosc/install/include \
		-DBLOSC_LIBRARY=$externalDependencies/c-blosc/install/lib/libblosc.so \
		-DCBENCH_ENABLE_SZ=ON \
		-DSZ_INCLUDE_PATH=$externalDependencies/SZ/sz/include \
		-DSZ_LIBRARY=$externalDependencies/SZ/install/lib/libSZ.so \
		-DZLIB_LIBRARY=$externalDependencies/SZ/install/lib/libzlib.so \
		-DZSTD_LIBRARY=$externalDependencies/SZ/install/lib/libzstd.so \
		-DCBENCH_ENABLE_ZFP=ON \
		-DZFP_INCLUDE_PATH=$externalDependencies/zfp/install/include \
		-DZFP_LIBRARY=$externalDependencies/zfp/install/lib64/libzfp.so \
		-DCBENCH_ENABLE_FPZIP=ON \
		-DFPZIP_INCLUDE_PATH=$externalDependencies/fpzip-1.2.0/inc/ \
		-DFPZIP_LIBRARY=$externalDependencies/fpzip-1.2.0/lib/libfpzip.a \
		-DCBENCH_ENABLE_ISABELA=ON \
		-DISABELA_INCLUDE_PATH=$externalDependencies/ISABELA-compress-0.2.1/include \
		-DISABELA_LIBRARY=$externalDependencies/ISABELA-compress-0.2.1/lib/libisabela.a \
		-DCBENCH_ENABLE_MGARD=ON \
		-DMGARD_INCLUDE_PATH=$externalDependencies/MGARD/install/include \
		-DMGARD_LIBRARY=$externalDependencies/MGARD/install/lib/libmgard.so \
		-DCBENCH_ENABLE_LOSSY_WAVE=ON \
		-DLOSSYWAVE_INCLUDE_PATH=$externalDependencies/VizAly-LossyWave/install/include \
		-DLOSSYWAVE_LIBRARY=$externalDependencies/VizAly-LossyWave/install/lib/liblossywave.so \
		-DLOSSYWAVE_LZ4_LIBRARY=$externalDependencies/VizAly-LossyWave/build/3rdparty/lz4/lz4-external/src/lz4-external/lib/liblz4.a \
		-DCBENCH_ENABLE_BIG_CRUNCH=ON \
		-DBIGCRUNCH_INCLUDE_PATH=$externalDependencies/VizAly-BigCrunch/install/include \
		-DBIGCRUNCH_LIBRARY=$externalDependencies/VizAly-BigCrunch/install/lib/libbigcrunch.so \
		-DCMAKE_BUILD_TYPE=$buildType
	fi
	
elif [ $buildOpt = "osx" ]; then
	echo "Building for osx ..."

	# OSX for Hoby
	cmake ../CBench \
		-DCMAKE_BUILD_TYPE=Debug \
		-DCBENCH_ENABLE_BLOSC=ON \
		-DBLOSC_INCLUDE_PATH=$externalDependencies/c-blosc/install/include \
		-DBLOSC_LIBRARY=$externalDependencies/c-blosc/install/lib/libblosc.a \
		-DCBENCH_ENABLE_SZ=ON \
		-DSZ_INCLUDE_PATH=$externalDependencies/SZ/sz/include \
		-DSZ_LIBRARY=$externalDependencies/SZ/install/lib/libSZ.dylib \
		-DZLIB_LIBRARY=$externalDependencies/SZ/install/lib/libzlib.dylib \
		-DZSTD_LIBRARY=$externalDependencies/SZ/install/lib/libzstd.dylib \
		-DCBENCH_ENABLE_BIG_CRUNCH=ON \
		-DBIGCRUNCH_INCLUDE_PATH=$externalDependencies/VizAly-BigCrunch/install/include \
		-DBIGCRUNCH_LIBRARY=$externalDependencies/VizAly-BigCrunch/install/lib/libbigcrunch.dylib \
		-DCBENCH_ENABLE_LOSSY_WAVE=ON \
		-DLOSSYWAVE_INCLUDE_PATH=$externalDependencies/VizAly-LossyWave/install/include \
		-DLOSSYWAVE_LIBRARY=$externalDependencies/VizAly-LossyWave/install/lib/liblossywave.dylib \
		-DLOSSYWAVE_LZ4_LIBRARY=$externalDependencies/VizAly-LossyWave/build/3rdparty/lz4/lz4-external/src/lz4-external/lib/liblz4.a \
		-DCBENCH_ENABLE_MGARD=OFF \
		-DMGARD_INCLUDE_PATH=$externalDependencies/MGARD/install/include \
		-DMGARD_LIBRARY=$externalDependencies/MGARD/install/lib/libmgard.dylib \
		-DCBENCH_ENABLE_ZFP=ON \
		-DZFP_INCLUDE_PATH=$externalDependencies/zfp/install/include \
		-DZFP_LIBRARY=$externalDependencies/zfp/install/lib/libzfp.dylib \
		-DCBENCH_ENABLE_NYX_LOADER=ON \
		-DHDF5_DIR=$externalDependencies/hdf5/install/share/cmake/hdf5 \
		-DCBENCH_ENABLE_FPZIP=ON \
		-DFPZIP_INCLUDE_PATH=$externalDependencies/fpzip-1.2.0/inc/ \
		-DFPZIP_LIBRARY=$externalDependencies/fpzip-1.2.0/lib/libfpzip.a \
		-DCBENCH_ENABLE_ISABELA=ON \
		-DISABELA_INCLUDE_PATH=$externalDependencies/ISABELA-compress-0.2.1/include \
		-DISABELA_LIBRARY=$externalDependencies/ISABELA-compress-0.2.1/lib/libisabela.a

elif [ $buildOpt = "gpu" ]; then
	echo "Building with gpu dependencies ..."

	cmake ../CBench $opt\
	    -DCBENCH_ENABLE_NYX_LOADER=ON \
		-DHDF5_DIR=$projectPath/ExternalDependencies/hdf5/install/share/cmake/hdf5 \
		-DCBENCH_ENABLE_SZ_GPU=ON \
		-DSZ_GPU_INCLUDE_PATH=$projectPath/ExternalDependencies/SZ-generic/sz/include \
		-DSZ_GPU_LIBRARY=$projectPath/ExternalDependencies/SZ-generic/install/lib/libSZ.so \
		-DZLIB_GPU_LIBRARY=$projectPath/ExternalDependencies/SZ-generic/install/lib/libzlib.so \
		-DZSTD_GPU_LIBRARY=$projectPath/ExternalDependencies/SZ-generic/install/lib/libzstd.so \
		-DCBENCH_ENABLE_ZFP_GPU=ON \
		-DZFP_GPU_INCLUDE_PATH=$projectPath/ExternalDependencies/gpu_zfp/install/include \
		-DZFP_GPU_LIBRARY=$projectPath/ExternalDependencies/gpu_zfp/install/lib64/libzfp.so \
		-DCBENCH_ENABLE_FPZIP=ON \
		-DFPZIP_INCLUDE_PATH=$projectPath/ExternalDependencies/fpzip-1.2.0/inc/ \
		-DFPZIP_LIBRARY=$projectPath/ExternalDependencies/fpzip-1.2.0/lib/libfpzip.a \
		-DCMAKE_BUILD_TYPE=$buildType

else
	echo "Build type " $buildOpt " not supported yet!"
fi

make -j


#
# Build Arguments
# --path <foder name> : build exec in that folder
# --externalDependencies <foder name>: path where external dependencies was built"

# --cooley : build for cooley
# --cori : build for cori at nersc

# --all : build with all options on
# --min : build with minimal options on
# --hacc : build for HACC only; no HDF5
# --gpu : build for gpu only

# --release: use release mode instead of debug

# --cleanAll <folder>: remove build folder and ExternalDependencies
# --clean <folder>: remove build folder

# -h or --help : help

# No argument is default:
# 	build path is build
# 	build type is Debug

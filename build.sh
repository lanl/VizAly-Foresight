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
		echo "  --externalDependencies <folder>: path where external dependencies was built"
		echo ""
		echo "  --gpu : makes a GPU build"
		echo ""
		echo "  --all : build with all options on"
		echo "  --min : build with minimal options on"
		echo "  --hacc: build for HACC only; no HDF5"
		echo ""
		echo "  --cori   : build for cori at nersc"
		echo ""
		echo "  --release: use release mode instead of debug"
		echo ""
		echo "  --clean <folder>: remove build folder"
		echo "  --cleanAll <folder>: remove build folder and ExternalDependencies"
		echo ""
		echo "  No arguments: CPU build, bild path is ./build, externalDependencies path is ./ExternalDependencies, type is Debug, "
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

	# clean
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

	# cleanAll
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

		echo "removing external dependencies"
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
pushd $buildDir 


# build
if [ $buildOpt = "min" ]; then
	echo "*** Default with minimal options ..."
	cmake ../CBench \
		-DCMAKE_BUILD_TYPE=$buildType

elif [ $buildOpt = "Default" ]; then
	echo "*** Default build on ..."
	if [ $buildPlatform = "cori" ]; then
		echo "... building on Cori ... not tested for a long while!!!"
		cmake ../CBench $opt\
			-DCBENCH_LOADER_GDA=ON \
			-DCBENCH_LOADER_NYX=ON \
			-DCBENCH_LOADER_GENERICBINARY=ON \
			-DHDF5_DIR=$externalDependencies/hdf5/install/share/cmake/hdf5 \
			-DCBENCH_COMPRESSOR_BLOSC=ON \
			-DBLOSC_INCLUDE_PATH=$externalDependencies/c-blosc/install/include \
			-DBLOSC_LIBRARY=$externalDependencies/c-blosc/install/lib/libblosc.so \
			-DCBENCH_COMPRESSOR_SZ=ON \
			-DSZ_INCLUDE_PATH=$externalDependencies/SZ/sz/include \
			-DSZ_LIBRARY=$externalDependencies/SZ/install/lib64/libSZ.a \
			-DZSTD_LIBRARY=$externalDependencies/SZ/install/lib64/libzstd.a \
			-DCBENCH_COMPRESSOR_ZFP=ON \
			-DZFP_INCLUDE_PATH=$externalDependencies/zfp/install/include \
			-DZFP_LIBRARY=$externalDependencies/zfp/install/lib64/libzfp.a \
			-DCBENCH_COMPRESSOR_FPZIP=ON \
			-DFPZIP_INCLUDE_PATH=$externalDependencies/fpzip-1.2.0/inc/ \
			-DFPZIP_LIBRARY=$externalDependencies/fpzip-1.2.0/lib/libfpzip.a \
			-DCMAKE_BUILD_TYPE=$buildType
	else
		export LD_LIBRARY_PATH=$externalDependencies/hdf5/install/lib:$LD_LIBRARY_PATH
		export PATH=$externalDependencies/hdf5/install/bin:$PATH

		cmake ../CBench $opt\
			-DCBENCH_LOADER_GDA=ON \
			-DCBENCH_LOADER_GENERICBINARY=ON \
			-DCBENCH_LOADER_NYX=ON \
			-DHDF5_DIR=$externalDependencies/hdf5/install/share/cmake/hdf5 \
			-DCBENCH_COMPRESSOR_BLOSC=ON \
			-DBLOSC_INCLUDE_PATH=$externalDependencies/c-blosc/install/include \
			-DBLOSC_LIBRARY=$externalDependencies/c-blosc/install/lib/libblosc.so \
			-DCBENCH_COMPRESSOR_ZFP=ON \
			-DZFP_INCLUDE_PATH=$externalDependencies/zfp/install/include \
			-DZFP_LIBRARY=$externalDependencies/zfp/install/lib64/libzfp.so \
			-DCBENCH_COMPRESSOR_FPZIP=ON \
			-DFPZIP_INCLUDE_PATH=$externalDependencies/fpzip-1.2.0/inc/ \
			-DFPZIP_LIBRARY=$externalDependencies/fpzip-1.2.0/lib/libfpzip.a \
			-DCBENCH_COMPRESSOR_SZ=ON \
			-DSZ_INCLUDE_PATH=$externalDependencies/SZ/sz/include \
			-DSZ_LIBRARY=$externalDependencies/SZ/install/lib64/libSZ.so \
			-DCBENCH_COMPRESSOR_MGARD=ON \
			-DMGARD_INCLUDE_PATH=$externalDependencies/MGARD/install/include \
			-DMGARD_LIBRARY=$externalDependencies/MGARD/install/lib64/libmgard.so \
			-DCMAKE_BUILD_TYPE=$buildType
	fi

elif [ $buildOpt = "hacc" ]; then
	echo "*** Building for HACC, no hdf5 ..."

	cmake ../CBench $opt\
		-DCBENCH_COMPRESSOR_BLOSC=ON \
		-DBLOSC_INCLUDE_PATH=$externalDependencies/c-blosc/install/include \
		-DBLOSC_LIBRARY=$externalDependencies/c-blosc/install/lib/libblosc.so \
		-DCBENCH_COMPRESSOR_SZ=ON \
		-DSZ_INCLUDE_PATH=$externalDependencies/SZ/sz/include \
		-DSZ_LIBRARY=$externalDependencies/SZ/install/lib64/libSZ.so \
		-DCBENCH_COMPRESSOR_ZFP=ON \
		-DZFP_INCLUDE_PATH=$externalDependencies/zfp/install/include \
		-DZFP_LIBRARY=$externalDependencies/zfp/install/lib64/libzfp.so \
		-DCBENCH_COMPRESSOR_FPZIP=ON \
		-DFPZIP_INCLUDE_PATH=$externalDependencies/fpzip-1.2.0/inc/ \
		-DFPZIP_LIBRARY=$externalDependencies/fpzip-1.2.0/lib/libfpzip.a \
		-DCMAKE_BUILD_TYPE=$buildType

elif [ $buildOpt = "all" ]; then
	echo "*** Building with all dependencies ..."

	if [ "$PLATFORM" = "travis" ]; then
	  echo "Travis: Using internal HDF5 build"
	  cmake ../CBench $opt\
		-DCBENCH_LOADER_GDA=ON \
		-DCBENCH_LOADER_GENERICBINARY=ON \
		-DCBENCH_COMPRESSOR_BLOSC=ON \
		-DBLOSC_INCLUDE_PATH=$externalDependencies/c-blosc/install/include \
		-DBLOSC_LIBRARY=$externalDependencies/c-blosc/install/lib/libblosc.so \
		-DCBENCH_COMPRESSOR_SZ=ON \
		-DSZ_INCLUDE_PATH=$externalDependencies/SZ/sz/include \
		-DSZ_LIBRARY=$externalDependencies/SZ/install/lib64/libSZ.so \
		-DCBENCH_COMPRESSOR_ZFP=ON \
		-DZFP_INCLUDE_PATH=$externalDependencies/zfp/install/include \
		-DZFP_LIBRARY=$externalDependencies/zfp/install/lib64/libzfp.so \
		-DCBENCH_COMPRESSOR_FPZIP=ON \
		-DFPZIP_INCLUDE_PATH=$externalDependencies/fpzip-1.2.0/inc/ \
		-DFPZIP_LIBRARY=$externalDependencies/fpzip-1.2.0/lib/libfpzip.a \
		-DCBENCH_COMPRESSOR_ISABELA=ON \
		-DISABELA_INCLUDE_PATH=$externalDependencies/ISABELA-compress-0.2.1/include \
		-DISABELA_LIBRARY=$externalDependencies/ISABELA-compress-0.2.1/lib/libisabela.a \
		-DCBENCH_COMPRESSOR_MGARD=ON \
		-DMGARD_INCLUDE_PATH=$externalDependencies/MGARD/install/include \
		-DMGARD_LIBRARY=$externalDependencies/MGARD/install/lib64/libmgard.so \
		-DCMAKE_BUILD_TYPE=$buildType
	else
	  cmake ../CBench $opt\
		-DCBENCH_LOADER_GDA=ON \
		-DCBENCH_LOADER_GENERICBINARY=ON \
		-DCBENCH_LOADER_NYX=ON \
		-DHDF5_DIR=$externalDependencies/hdf5/install/share/cmake/hdf5 \
		-DCBENCH_COMPRESSOR_BLOSC=ON \
		-DBLOSC_INCLUDE_PATH=$externalDependencies/c-blosc/install/include \
		-DBLOSC_LIBRARY=$externalDependencies/c-blosc/install/lib/libblosc.so \
		-DCBENCH_COMPRESSOR_SZ=ON \
		-DSZ_INCLUDE_PATH=$externalDependencies/SZ/sz/include \
		-DSZ_LIBRARY=$externalDependencies/SZ/install/lib64/libSZ.so \
		-DCBENCH_COMPRESSOR_ZFP=ON \
		-DZFP_INCLUDE_PATH=$externalDependencies/zfp/install/include \
		-DZFP_LIBRARY=$externalDependencies/zfp/install/lib64/libzfp.so \
		-DCBENCH_COMPRESSOR_FPZIP=ON \
		-DFPZIP_INCLUDE_PATH=$externalDependencies/fpzip-1.2.0/inc/ \
		-DFPZIP_LIBRARY=$externalDependencies/fpzip-1.2.0/lib/libfpzip.a \
		-DCBENCH_COMPRESSOR_ISABELA=ON \
		-DISABELA_INCLUDE_PATH=$externalDependencies/ISABELA-compress-0.2.1/include \
		-DISABELA_LIBRARY=$externalDependencies/ISABELA-compress-0.2.1/lib/libisabela.a \
		-DCBENCH_COMPRESSOR_MGARD=ON \
		-DMGARD_INCLUDE_PATH=$externalDependencies/MGARD/install/include \
		-DMGARD_LIBRARY=$externalDependencies/MGARD/install/lib64/libmgard.so \
		-DCMAKE_BUILD_TYPE=$buildType
	fi

elif [ $buildOpt = "gpu" ]; then
	echo "*** Building with gpu dependencies ..."

	cmake ../CBench $opt\
		-DCBENCH_LOADER_GDA=ON \
		-DCBENCH_LOADER_NYX=ON \
		-DHDF5_DIR=$projectPath/ExternalDependencies/hdf5/install/share/cmake/hdf5 \
		-DCBENCH_COMPRESSOR_SZ_GPU=ON \
		-DSZ_GPU_INCLUDE_PATH=$projectPath/ExternalDependencies/SZ-generic/sz/include \
		-DSZ_GPU_LIBRARY=$projectPath/E
		-DZFP_GPU_INCLUDE_PATH=$projectPath/ExternalDependencies/gpu_zfp/install/include \
		-DZFP_GPU_LIBRARY=$projectPath/ExternalDependencies/gpu_zfp/install/lib64/libzfp.so \
		-DCBENCH_COMPRESSOR_FPZIP=ON \
		-DFPZIP_INCLUDE_PATH=$projectPath/ExternalDependencies/fpzip-1.2.0/inc/ \
		-DFPZIP_LIBRARY=$projectPath/ExternalDependencies/fpzip-1.2.0/lib/libfpzip.a \
		-DCMAKE_BUILD_TYPE=$buildType

else
	echo "*** Build type " $buildOpt " not supported yet!"
fi

make -j

popd


#
# Build Arguments
# --path <foder name> : build exec in that folder
# sh  <foder name>: path where external dependencies was built"

# --cori : special build for cori at nersc

# --all : build with all options on
# --min : build with minimal options on
# --hacc : build for HACC only; no HDF5
# --gpu : build for gpu only

# --release: use release mode instead of debug

# --cleanAll <folder>: remove build folder and ExternalDependencies
# --clean <folder>: remove build folder

# -h or --help : help

# No argument is default:
#   will build on CPU
# 	build path is build
#   external dependencies path is ExternalDependencies
# 	build type is Debug

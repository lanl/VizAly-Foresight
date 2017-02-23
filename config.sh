# source environment
source env/bash.darwin.gcc-openmp

# input variables
build_type=${1:-Release}

# move into the 3rdparty folder
root=$(eval pwd)
cd source/3rdparty

# configure sz
echo " "
echo ===============================================================================================
echo starting sz configuration
echo ===============================================================================================
echo " "
cd sz-1.4.9-hacc
if [ "$build_type" == "Release" ]
then
	./configure --prefix=$(eval pwd)/build
else
	./configure --prefix=$(eval pwd)/build --enable-debug
fi
cd ..

# configure xz/lzma
echo " "
echo ===============================================================================================
echo starting xz/lzma configuration
echo ===============================================================================================
echo " "
cd xz-5.2.3
if [ "$build_type" == "Release" ]
then
	./configure --prefix=$(eval pwd)/build
else
	./configure --prefix=$(eval pwd)/build --enable-debug
fi
cd ..

# run cmake
echo " "
echo ===============================================================================================
echo starting cmake with build type $build_type
echo ===============================================================================================
echo " "
cd $root
mkdir -p binary
cd binary
cmake ../source -DCMAKE_INSTALL_PREFIX=../build -DBUILD_BENCHMARKS=OFF -DBUILD_STATIC=OFF -DBUILD_TESTS=OFF -DCMAKE_BUILD_TYPE=$build_type -DDATA_LOADER_BUILD_STATIC=ON -DDATA_LOADER_INSTALL=OFF -DGENERICIO_BUILD_STATIC=ON -DUTILITY_BUILD_STATIC=ON -DUTILITY_INSTALL=OFF
cd $root

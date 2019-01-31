#!/bin/bash

# MGARD
echo "Building MGARD ... "

git clone https://github.com/CODARcode/MGARD.git
cd MGARD
# Jan 2019 hash
git checkout a6db5d03e59c156b1e2df32db845a4c38e74b04a

cd include
if [ ! -f mgard_capi.h.bak ]; then
	echo "Patching mgard_capi.h"
	sed -i.bak "40a\\extern unsigned char \*mgard\_compress\(int itype\_flag, void  \*data, int \&out\_size, int nrow, int ncol, int nfib, void\* tol_in\);" mgard_capi.h
	sed -i "39d" mgard_capi.h
fi
cd ..

cd src
if [ ! -f mgard_capi.cpp.bak ]; then
	echo "Patching mgard_capi.cpp"
	sed -i.bak "s/\"C\"//" mgard_capi.cpp
fi
cd ..

mkdir install
mkdir install/lib
mkdir install/include
mkdir build

if [ ! -f CMakeLists.txt.bak ]; then
	echo "Patching CMakeLists.txt"
	sed -i.bak "s|add\\_library[(][^$]*|&SHARED\\ |g" CMakeLists.txt
fi

cd build

g++ -std=c++11 -Wall -shared -fPIC -o libmgard.so ../src/mgard_capi.cpp ../src/mgard_nuni.cpp ../src/mgard.cpp -I../include

echo "Installing MGARD"
cp libmgard.so ../install/lib
cp ../include/* ../install/include/
#cmake .. -DCMAKE_INSTALL_PREFIX=../install
#make -j
#make install
cd ..
cd ..

echo "Building MGARD done!"

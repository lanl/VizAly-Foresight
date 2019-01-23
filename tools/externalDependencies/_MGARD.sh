#!/bin/bash

# MGARD
echo "Building MGARD ... "

git clone https://github.com/CODARcode/MGARD.git
cd MGARD
# Jan 2019 hash
git checkout a6db5d03e59c156b1e2df32db845a4c38e74b04a
mkdir install
mkdir build

if [ ! -f CMakeLists.txt.bak ]; then
	echo "Patching CMakeLists.txt"
	sed -i.bak "s|add\\_library[(][^$]*|&SHARED\\ |g" CMakeLists.txt
fi

cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j
make install
cd ..
cd ..

echo "Building MGARD done!"

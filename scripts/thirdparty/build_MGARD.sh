#!/bin/bash

#!/bin/bash

# MGARD
echo "Building MGARD ... "

git clone https://github.com/CODARcode/MGARD.git
cd MGARD
git checkout 03c3e6e01e3654b11af5d7d4fbeebd5ae59ca84f	# This one works

mkdir install
mkdir build

cd build

echo "Building MGARD"
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make -j

echo "Installing MGARD"
make install
cd ..
cd ..

echo "Building MGARD done!"
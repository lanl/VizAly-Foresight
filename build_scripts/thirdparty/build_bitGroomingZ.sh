#!/bin/bash

echo "Building BitGroomingZ ... "

git clone https://github.com/disheng222/BitGroomingZ.git
cd BitGroomingZ/

mkdir install

./configure --prefix=$(pwd)/install
make -j
make install

cd ..

echo "Building BitGroomingZ master done!"
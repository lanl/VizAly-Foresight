#!/bin/bash

# ISABELA
echo "Building ISABELA ... "

wget http://freescience.org/cs/ISABELA/code/ISABELA-compress-0.2.1.tar.gz
tar -zxvf ISABELA-compress-0.2.1.tar.gz
cd ISABELA-compress-0.2.1
make -j
cd ..

echo "Building ISABELA done!"


#!/bin/bash

# ISABELA
echo "Building ISABELA ... "

echo "Needs libssl-dev for md5"
echo "Needs gsl for wavelets"

wget http://freescience.org/cs/ISABELA/code/ISABELA-compress-0.2.1.tar.gz
tar -zxvf ISABELA-compress-0.2.1.tar.gz
cd ISABELA-compress-0.2.1
make -j

if [ $? == 1 ]; then
echo "Building ISABELA done!"
else
echo "Building ISABELA failed!"
fi

cd ..

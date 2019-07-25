#!/bin/bash

# ISABELA
echo "Building ISABELA ... "

#Needs libssl-dev for md5
#Needs gsl for wavelets

wget http://freescience.org/cs/ISABELA/code/ISABELA-compress-0.2.1.tar.gz
tar -zxvf ISABELA-compress-0.2.1.tar.gz
cd ISABELA-compress-0.2.1
make -j
cd ..

echo "Building ISABELA done!"


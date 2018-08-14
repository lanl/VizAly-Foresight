#!/bin/bash

rootDir=`pwd`

git pull

if [ ! -d Z-checker ]; then
	echo "Error: no Z-checker directory."
	echo "Please use z-checker-installer.sh to perform the installation first."
	exit
fi

#---------- download Z-checker --------------
cd Z-checker
git pull
make
make install
cp ../zc-patches/generateReport.sh ./examples/

cd examples
make clean
make

#---------- download ZFP and set the configuration -----------
cd $rootDir
cd zfp
git pull

cd -
cp zfp-patches/zfp-zc.c zfp/utils
cp zfp-patches/*.sh zfp/utils

make

#---------- download SZ and set the configuration -----------
cd $rootDir
cd SZ
git pull
make
make install

cd example
cp ../../Z-checker/examples/zc.config .
cp ../../sz-patches/sz-zc-ratedistortion.sh .
cp ../../sz-patches/testfloat_CompDecomp.sh .
cp ../../sz-patches/testdouble_CompDecomp.sh .

make clean
make

#!/bin/bash

if [ $# != 1 ]
then
        echo Usage: please specify a unique case directory name.
        echo Example: $0 case1 
        exit
fi

caseName=$1

rootDir=`pwd`

if [ ! -d Z-checker ]; then
	echo "Error: missing Z-checker directory. "
	echo "Please run z-checker-install.sh first."
	exit
fi

echo Create a new case $caseName for Z-checker
if [ -d Z-checker/${caseName}-pwr ]; then
	echo "Conflict: the case ${caseName}-pwr already exists."
	echo "Please remove the existing case using removeZCCase.sh before creating another one with the same case name."
	exit
fi
cd Z-checker
./createNewCase.sh ${caseName}-pwr

echo "Create a new case (fast mode) for SZ"
cd ../SZ
sz_caseName=${caseName}-pwr_fast
if [ ! -d $sz_caseName ]; then
	mkdir $sz_caseName
fi
cp example/sz-zc-ratedistortion.sh $sz_caseName
cp example/testfloat_CompDecomp.sh $sz_caseName
cp example/testdouble_CompDecomp.sh $sz_caseName
cp example/zc.config $sz_caseName
cp ../sz-patches/sz.config.fast_mode $sz_caseName/sz.config
cd $sz_caseName
ln -s "$rootDir/SZ/example/testfloat_CompDecomp" testfloat_CompDecomp
patch -p0 < ../../sz-patches/testfloat_CompDecomp_fast.sh.patch
ln -s "$rootDir/SZ/example/testdouble_CompDecomp" testdouble_CompDecomp
patch -p0 < ../../sz-patches/testdouble_CompDecomp_fast.sh.patch
cd ..

echo "Create a new case (default mode) for SZ"
cd ../SZ
sz_caseName=${caseName}-pwr_deft
if [ ! -d $sz_caseName ]; then
	mkdir $sz_caseName
fi
cp example/sz-zc-ratedistortion.sh $sz_caseName
cp example/testfloat_CompDecomp.sh $sz_caseName
cp example/testdouble_CompDecomp.sh $sz_caseName
cp example/zc.config $sz_caseName
cp ../sz-patches/sz.config.default_mode $sz_caseName/sz.config
cd $sz_caseName
ln -s "$rootDir/SZ/example/testfloat_CompDecomp" testfloat_CompDecomp
patch -p0 < ../../sz-patches/testfloat_CompDecomp_deft.sh.patch
ln -s "$rootDir/SZ/example/testdouble_CompDecomp" testdouble_CompDecomp
patch -p0 < ../../sz-patches/testdouble_CompDecomp_deft.sh.patch
cd ..

echo Create a new case for ZFP
cd ../zfp
zfp_caseName=${caseName}-p
if [ ! -d $zfp_caseName ]; then
	mkdir $zfp_caseName
fi
cp utils/*.sh $zfp_caseName

echo Modify Z-checker/$caseName-pwr/zc.config
cd ../Z-checker/$caseName-pwr
./modifyZCConfig zc.config compressors "sz_f:../../SZ/${caseName}-pwr_fast sz_d:../../SZ/${caseName}-pwr_deft zfp:../../zfp/${zfp_caseName}"

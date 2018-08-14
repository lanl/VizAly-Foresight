#!/bin/bash

if [[ $# < 3 ]]; then
	echo "Usage: $0 [datatype (-f or -d)] [errorBoundMode] [data directory] [dimension sizes....]"
	echo Example: $0 -f ABS /home/shdi/CESM-testdata/1800x3600 3600 1800
	exit
fi

datatype=$1
errBoundMode=$2
dataDir="$3"
dim1=$4
dim2=$5
dim3=$6
dim4=$7

#Note: If you run this script by z-checker-installer, ZFP_Err_Bounds will be overwritten by ../../errBounds.cfg as follows.
ZFP_Err_Bounds="1E-1 1E-2 1E-3 1E-4"

echo $errBoundMode
if [ -f ../../errBounds.cfg ]; then
	if [[ $errBoundMode == "PW_REL" ]]; then
		zfp_err_env="`cat ../../errBounds_pwr.cfg | grep -v "#" | grep ZFP_ERR_BOUNDS`"
	else
		zfp_err_env="`cat ../../errBounds.cfg | grep -v "#" | grep ZFP_ERR_BOUNDS`"
	fi
	echo "export $zfp_err_env" > env.tmp
	source env.tmp
	rm env.tmp
	ZFP_Err_Bounds="`echo $ZFP_ERR_BOUNDS`"
	echo $ZFP_Err_Bounds
fi

for errBound in $ZFP_Err_Bounds
do
	./zfp-zc-dir.sh $datatype $errBoundMode $errBound "$dataDir" $dim1 $dim2 $dim3 $dim4
done


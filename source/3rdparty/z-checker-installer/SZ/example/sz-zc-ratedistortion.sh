#!/bin/bash

if [[ $# < 4 ]]; then
	echo "Usage: $0 [data type (-f or -d)] [errBoundMode] [data directory] [dimension sizes....]"
	echo "Example: $0 -f ABS /home/shdi/CESM-testdata/1800x3600 3600 1800"
	exit
fi

datatype=$1
errBoundMode=$2
dataDir=$3
#dataDir=/home/fti/SZ_C_version/CESM-testdata/1800x3600

dim1=$4
dim2=$5
dim3=$6
dim4=$7

#Note: If you run this script by z-checker-installer, SZ_Err_Bounds will be overwritten by ../../errBounds.cfg as follows.
SZ_Err_Bounds="1E-1 1E-2 1E-3 1E-4"

if [ -f ../../errBounds.cfg ]; then
	if [[ $errBoundMode == "PW_REL" ]];then
		sz_err_env="`cat ../../errBounds_pwr.cfg | grep -v "#" | grep SZ_ERR_BOUNDS`"
	else
		sz_err_env="`cat ../../errBounds.cfg | grep -v "#" | grep SZ_ERR_BOUNDS`"
	fi
	echo "export $sz_err_env" > env.tmp
	source env.tmp
	rm env.tmp
	SZ_Err_Bounds="`echo $SZ_ERR_BOUNDS`"
fi

for errBound in $SZ_Err_Bounds
do
	if [[ $datatype == "-f" ]]; then
		./testfloat_CompDecomp.sh $errBoundMode $errBound "$dataDir" $dim1 $dim2 $dim3 $dim4
	elif [[ $datatype == "-d" ]]; then
		./testdouble_CompDecomp.sh $errBoundMode $errBound "$dataDir" $dim1 $dim2 $dim3 $dim4
	else
		echo "Error: datatype = $datatype . "
		echo "Note: datatype can only be either -f or -d."
		exit
	fi
done


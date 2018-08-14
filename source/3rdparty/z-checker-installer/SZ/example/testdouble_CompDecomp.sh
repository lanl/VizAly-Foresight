#!/bin/bash

if [[ $# < 4 ]]
then
	echo Usage: $0 [errBoundMode] [error bound] [data directory] [dimension sizes....]
	echo Example: $0 ABS 1E-4 /home/fti/SZ_C_version/CESM-testdata/1800x3600 3600 1800
	exit
fi

errBoundMode=$1
absErrBound=$2
dataDir=$3
dim1=$4
dim2=$5
dim3=$6
dim4=$7

compressor=sz

fileList=`cd "$dataDir";ls *.dat`
for file in $fileList
do
        echo testdouble_CompDecomp sz.config zc.config "${compressor}($absErrBound)" "$file" $errBoundMode $absErrBound "$dataDir/$file" $dim1 $dim2 $dim3 $dim4
        ./testdouble_CompDecomp sz.config zc.config "${compressor}($absErrBound)" "$file" $errBoundMode $absErrBound "$dataDir/$file" $dim1 $dim2 $dim3 $dim4
done

echo "complete"


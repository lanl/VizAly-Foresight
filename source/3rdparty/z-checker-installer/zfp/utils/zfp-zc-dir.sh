#!/bin/bash

if [[ $# < 4 ]]
then
	echo "Usage: $0 [datatype (-f or -d)] [errBoundMode] [error bound] [data directory] [dimension sizes....]"
	echo "Example: $0 -f ABS 1E-4 /home/fti/SZ_C_version/CESM-testdata/1800x3600 3600 1800"
	exit
fi

datatype=$1
errBoundMode=$2
absErrBound=$3
dataDir="$4"
dim1=$5
dim2=$6
dim3=$7
dim4=$8

fileList=`cd "$dataDir";ls *.dat`
for file in $fileList
do
	echo ./zfp-zc.sh $datatype $errBoundMode ${dataDir}/${file} ${file} $absErrBound $dim1 $dim2 $dim3 $dim4
	./zfp-zc.sh $datatype $errBoundMode "${dataDir}/${file}" ${file} $absErrBound $dim1 $dim2 $dim3 $dim4
done

echo "complete"


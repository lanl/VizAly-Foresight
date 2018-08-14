#!/bin/bash

if [[ $# < 5 ]]
then
	echo "Usage: $0 [datatype (-f or -d)] [errBoundMode] [dataFilePath] [varName] [errBound] [dimension sizes....]"
	echo Example: $0 -f [ABS/REL/PW_REL] CESM-testdata/CLDLOW_1_1800_3600.dat CLDLOW 1E-4 3600 1800
	exit
fi

cmdDir=../bin

datatype=$1
errBoundMode=$2
dataFilePath=$3
varName=$4
errBound=$5
let dim=$#-5

if [[ $errBoundMode == "ABS" ]]; then
	if [[ $dim == 1 ]]
	then
		echo ${cmdDir}/zfp-zc -s $datatype -a ${errBound} -${dim} $6 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}"
		${cmdDir}/zfp-zc -s $datatype -a ${errBound} -${dim} $6 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}"
	elif [[ $dim == 2 ]]
	then
		echo ${cmdDir}/zfp-zc -s $datatype -a ${errBound} -${dim} $6 $7 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}"
		${cmdDir}/zfp-zc -s $datatype -a ${errBound} -${dim} $6 $7 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}"
	elif [[ $dim == 3 ]]
	then
		echo ${cmdDir}/zfp-zc -s $datatype -a ${errBound} -${dim} $6 $7 $8 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}"
		${cmdDir}/zfp-zc -s $datatype -a ${errBound} -${dim} $6 $7 $8 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}"
	fi
elif [[ $errBoundMode == "REL" ]]; then
	if [[ $dim == 1 ]]
	then
		echo ${cmdDir}/zfp-zc -s $datatype -a ${errBound} -${dim} $6 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}" -l
		${cmdDir}/zfp-zc -s $datatype -a ${errBound} -${dim} $6 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}" -l
	elif [[ $dim == 2 ]]
	then
		echo ${cmdDir}/zfp-zc -s $datatype -a ${errBound} -${dim} $6 $7 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}" -l
		${cmdDir}/zfp-zc -s $datatype -a ${errBound} -${dim} $6 $7 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}" -l
	elif [[ $dim == 3 ]]
	then
		echo ${cmdDir}/zfp-zc -s $datatype -a ${errBound} -${dim} $6 $7 $8 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}" -l
		${cmdDir}/zfp-zc -s $datatype -a ${errBound} -${dim} $6 $7 $8 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}" -l
	fi
elif [[ $errBoundMode == "PW_REL" ]]; then
	if [[ $dim == 1 ]]
	then
		echo ${cmdDir}/zfp-zc -s $datatype -p ${errBound} -${dim} $6 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}"
		${cmdDir}/zfp-zc -s $datatype -p ${errBound} -${dim} $6 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}"
	elif [[ $dim == 2 ]]
	then
		echo ${cmdDir}/zfp-zc -s $datatype -p ${errBound} -${dim} $6 $7 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}"
		${cmdDir}/zfp-zc -s $datatype -p ${errBound} -${dim} $6 $7 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}"
	elif [[ $dim == 3 ]]
	then
		echo ${cmdDir}/zfp-zc -s $datatype -p ${errBound} -${dim} $6 $7 $8 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}"
		${cmdDir}/zfp-zc -s $datatype -p ${errBound} -${dim} $6 $7 $8 -i ${dataFilePath} -k "zfp(${errBound})" -v "${varName}"
	fi
fi

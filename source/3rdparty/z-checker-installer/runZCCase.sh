#!/bin/bash

datatype=$1
if [[ $# < 4 || ( $datatype != "-f" && $datatype != "-d" ) ]]
then
	echo "Usage: $0 [datatype (-f or -d)] [errBoundMode] [testcase] [data dir] [dimensions....]"
	echo "Example: $0 -f ABS testcase1 CESM-testdata/1800x3600 3600 1800"
	exit
fi 

errBoundMode=$2
testcase=$3
dataDir=`cd "$4"; pwd`
dim1=$5
dim2=$6
dim3=$7
dim4=$8

rootDir=`pwd`

if [[ $errBoundMode == "ABS" ]]; then
	if [ ! -d "Z-checker/$testcase" ]; then
		echo "Error: Testcase $testcase doesn't exist!"
		exit
	fi
elif [[ $errBoundMode == "REL" ]]; then
	if [ ! -d "Z-checker/$testcase" ]; then
		echo "Error: Testcase $testcase doesn't exist!"
		exit
	fi
elif [[ $errBoundMode == "PW_REL" ]]; then
	if [ ! -d "Z-checker/$testcase-pwr" ]; then
		echo "Error: Testcase $testcase for PW_REL doesn't exist!"
		exit
	fi
fi

if [ ! -d "$dataDir" ]; then
	echo "Error: $dataDir doesn't exist!"
	exit
fi

envConfigPath="$rootDir/Z-checker/examples/env_config.sh"
GNUPLOT_EXE_PATH=`which gnuplot`
if [ ! -x "$GNUPLOT_EXE_PATH" ]; then
	if [ -f $envConfigPath ]; then
		source $envConfigPath
	else
		echo "Error: gnuplot is not executable and cannot find Z-checker/examples/env_config.sh either."
		exit
	fi
fi

LATEXMK_EXE_PATH=`which latexmk`
if [ ! -x "$LATEXMK_EXE_PATH" ]; then
	if [ -z "$GNUPLOT_PATH" ]; then
		if [ -f $envConfigPath ]; then
			source $envConfigPath
		else
			echo "Error: latexmk is not executable and cannot find Z-checker/examples/env_config.sh either."
			exit
		fi
	fi
fi

if [[ $errBoundMode == "PW_REL" ]]; then
	cd SZ/${testcase}-pwr_fast
else
	cd SZ/${testcase}_fast
fi
echo ./sz-zc-ratedistortion.sh $datatype $errBoundMode $dataDir $dim1 $dim2 $dim3 $dim4
./sz-zc-ratedistortion.sh $datatype $errBoundMode $dataDir $dim1 $dim2 $dim3 $dim4

cd $rootDir
if [[ $errBoundMode == "PW_REL" ]]; then
	cd SZ/${testcase}-pwr_deft
else
	cd SZ/${testcase}_deft
fi
echo ./sz-zc-ratedistortion.sh $datatype $errBoundMode $dataDir $dim1 $dim2 $dim3 $dim4
./sz-zc-ratedistortion.sh $datatype $errBoundMode $dataDir $dim1 $dim2 $dim3 $dim4

cd $rootDir
if [[ $errBoundMode == "PW_REL" ]]; then
	cd zfp/${testcase}-p
else
	cd zfp/${testcase}
fi
echo ./zfp-zc-ratedistortion.sh $datatype $errBoundMode $dataDir $dim1 $dim2 $dim3 $dim4
./zfp-zc-ratedistortion.sh $datatype $errBoundMode $dataDir $dim1 $dim2 $dim3 $dim4

cd $rootDir
if [[ $errBoundMode == "PW_REL" ]]; then
	cd Z-checker/${testcase}-pwr
else
	cd Z-checker/${testcase}
fi
echo ./analyzeDataProperty.sh $datatype $dataDir $dim1 $dim2 $dim3 $dim4
./analyzeDataProperty.sh $datatype $dataDir $dim1 $dim2 $dim3 $dim4

if [[ $errBoundMode == "PW_REL" ]]; then
	sz_err_env="`cat ../../errBounds_pwr.cfg | grep -v "#" | grep comparisonCases`"
else
	sz_err_env="`cat ../../errBounds.cfg | grep -v "#" | grep comparisonCases`"
fi
echo "export $sz_err_env" > env.tmp
source env.tmp
rm env.tmp
SZ_Err_Bounds="`echo $comparisonCases`"

echo comparisonCases=$comparisonCases
./modifyZCConfig zc.config comparisonCases "$comparisonCases"

if [[ $errBoundMode == "PW_REL" ]]; then
	echo ./generateReport.sh ${testcase}-pwr
	./generateReport.sh "${testcase} with PW_REL"
else
	echo ./generateReport.sh ${testcase}
	./generateReport.sh "${testcase}"
fi

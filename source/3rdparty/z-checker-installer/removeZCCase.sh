#!/bin/bash

if [[ $# < 2 ]]; then
	echo Usage: $0 [errBoundMode] [testcase name]
	exit
fi

errBoundMode=$1
testcase=$2
if [[ $errBoundMode == "ABS" ]]; then
	if [ -d Z-checker/$testcase ]; then
		rm -rf Z-checker/$testcase
		rm -rf SZ/${testcase}_fast
		rm -rf SZ/${testcase}_deft
		rm -rf zfp/${testcase}
	else
		echo No such testcase: $testcase
		exit
	fi
elif [[ $errBoundMode == "PW_REL" ]]; then
	if [ -d Z-checker/${testcase}-pwr ]; then
		rm -rf Z-checker/$testcase-pwr
		rm -rf SZ/${testcase}-pwr_fast
		rm -rf SZ/${testcase}-pwr_def
		rm -rf zfp/${testcase}-pwr
	else
		echo No such testcase: $testcase-pwr
		exit
	fi
else
	echo Error: wrong errBoundMode.
	exit
fi

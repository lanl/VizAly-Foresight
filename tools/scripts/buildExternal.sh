#!/bin/bash

# Builds dependencies
currentDir=$(pwd)

for file in build_*.sh; 
do
	source $file; 
	cd $currentDir # to guard against build failure
done


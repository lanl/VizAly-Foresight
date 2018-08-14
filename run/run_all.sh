#!/bin/bash

# parameters
timestep=42
reps=10
ntasks=256

# evaluate command line parameters
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		-t|--timestep)
		timestep="$2"
		shift
		;;
		-r|--reps)
		reps="$2"
		shift
		;;
		-n|--ntasks)
		ntasks="$2"
		shift
		;;
	esac
	shift
done

# call job scripts
./run_blosc.sh --timestep $timestep --reps $reps --ntasks $ntasks
./run_lzma.sh --timestep $timestep --reps $reps --ntasks $ntasks
./run_fpzip.sh --timestep $timestep --reps $reps --ntasks $ntasks
./run_zfp.sh --timestep $timestep --reps $reps --ntasks $ntasks

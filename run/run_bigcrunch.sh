#!/bin/bash

# parameters
timestep=399
reps=10
ntasks=256

# evaluate command line parameters
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		--timestep)
		timestep="$2"
		shift
		;;
		--reps)
		reps="$2"
		shift
		;;
		--ntasks)
		ntasks="$2"
		shift
		;;
	esac
	shift
done

# init variables
#num_threads=( 1 2 4 8 16 )
num_threads=( 1 2 4 8 )
#in_name="/path/to/data/m000.full.mpicosmo.$timestep"

for nthreads in "${num_threads[@]}"
do
	for error in {-6..-1}
	do
		sbatch -n $ntasks -c $nthreads --ntasks-per-node="$(expr 72 / $nthreads)" bigcrunch_benchmark.sh --in_name $in_name --reps $reps --nthreads $nthreads --tolerance 0 --error $error
	done
done

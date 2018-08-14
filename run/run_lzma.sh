#!/bin/bash

# parameters
timestep=0
format="HACC"
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
		--format)
		format="$2"
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
presets=( 1 1e 2 2e 3 3e 4 4e 5 5e 6 6e 7 7e 8 8e 9 9e )
num_threads=( 1 2 4 8 16 )
#in_name="/path/to/data/m000.full.mpicosmo.$timestep"

for nthreads in "${num_threads[@]}"
do
	for preset in "${presets[@]}"
	do
		sbatch -n $ntasks -c $nthreads --ntasks-per-node="$(expr 72 / $nthreads)" lzma_benchmark.sh --in_name $in_name --format $format --reps $reps --preset $preset --nthreads $nthreads
	done
done

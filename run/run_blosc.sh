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
cnames=( blosclz zlib snappy lz4 lz4hc zstd )
shuffles=( NO_SHUFFLE SHUFFLE BIT_SHUFFLE )
num_threads=( 1 2 4 8 16 )
#in_name="/path/to/data/m000.full.mpicosmo.$timestep"

for nthreads in "${num_threads[@]}"
do
	for cname in "${cnames[@]}"
	do
		for shuffle in "${shuffles[@]}"
		do
			for clevel in {1..9}
			do
				echo sbatch -n $ntasks -c $nthreads --ntasks-per-node="$(expr 72 / $nthreads)" blosc_benchmark.sh --in_name $in_name --format $format --reps $reps --cname $cname --shuffle $shuffle --clevel $clevel --nthreads $nthreads
			done
		done
	done
done

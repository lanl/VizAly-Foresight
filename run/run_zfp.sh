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
num_threads=( 1 2 4 8 16 )
#in_name="/path/to/data/m000.full.mpicosmo.$timestep"

for nthreads in "${num_threads[@]}"
do
	for param in {0..7}
	do
		sbatch -n $ntasks -c $nthreads --ntasks-per-node="$(expr 72 / $nthreads)" zfp_benchmark.sh --in_name $in_name --format $format --reps $reps --mode ACCURACY --param $param --nthreads $nthreads
	done

	for param in {10..32}
	do
		sbatch -n $ntasks -c $nthreads --ntasks-per-node="$(expr 72 / $nthreads)" zfp_benchmark.sh --in_name $in_name --format $format --reps $reps --mode PRECISION --param $param --nthreads $nthreads
	done

	for param in {10..32}
	do
		sbatch -n $ntasks -c $nthreads --ntasks-per-node="$(expr 72 / $nthreads)" zfp_benchmark.sh --in_name $in_name --format $format --reps $reps --mode RATE --param $param --nthreads $nthreads
	done
done

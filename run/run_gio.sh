#!/bin/bash

# parameters
in_name=""
format=""
reps=10
ntasks=256

# evaluate command line parameters
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		--in_name)
		in_name="$2"
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

for nthreads in "${num_threads[@]}"
do
	for i in {1..16}
	do
		for j in {1..32}
		do
			num_buckets=$(echo "2^i" | bc -l)
			blocksize=$(echo "2^j" | bc -l)
			sbatch -n $ntasks -c $nthreads gio_benchmark.sh --in_name $in_name --format $format --reps $reps --num_buckets $num_buckets --blocksize $blocksize
		done
	done
done

#!/bin/bash

# parameters
timestep=499
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
sz_mode=( SZ_BEST_SPEED SZ_DEFAULT_COMPRESSION SZ_BEST_COMPRESSION )
gzip_mode=( Gzip_BEST_SPEED Gzip_BEST_COMPRESSION Gzip_DEFAULT_COMPRESSION )
#error_mode=( ABS_AND_REL ABS_OR_REL ABS REL PW_REL )
error_mode=( REL )

for nthreads in "${num_threads[@]}"
do
	for layers in {1..3}
	do
		for sm in "${sz_mode[@]}"
		do
			for em in "${error_mode[@]}"
			do
				for rel in {1..6}
				do
					if [ "$sm" == "SZ_BEST_COMPRESSION" ]
					then
						for gm in "${gzip_mode[@]}"
						do
							sbatch -n $ntasks -c $nthreads --ntasks-per-node="$(expr 72 / $nthreads)" sz_benchmark.sh --in_name $in_name --format $format --reps $reps --nthreads $nthreads --layers $layers --sz_mode $sm --error_bound_mode $em --rel_bound_ratio $(echo "10^-$rel" | bc -l) --gzip_mode $gm
						done
					else
						sbatch -n $ntasks -c $nthreads --ntasks-per-node="$(expr 72 / $nthreads)" sz_benchmark.sh --in_name $in_name --format $format --reps $reps --nthreads $nthreads --layers $layers --sz_mode $sm --error_bound_mode $em --rel_bound_ratio $(echo "10^-$rel" | bc -l)
					fi
				done
			done
		done
	done
done

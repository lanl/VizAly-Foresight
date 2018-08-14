#!/bin/bash
#SBATCH --output logs/blosc_benchmark.%j.out
#SBATCH --error logs/blosc_benchmark.%j.err
#SBATCH --partition=scaling

# init parameters
in_name=""
result_name=""
format=""
reps=10
cname="blosclz"
clevel=6
shuffle=1
nthreads=1

# evaluate command line parameters
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		-h|--help)
		echo "USAGE: sbatch blosc_benchmark.sh"
		echo " "
		echo "	-h | --help:"
		echo "		Shows this dialog."
		echo " "
		echo "	--in_name:"
		echo "		Input path to data."
		echo " "
		echo "	--format:"
		echo "		Data format (HACC, MPAS, EXALT)."
		echo " "
		echo "	--reps:"
		echo "		Number of repetitions for the benchmark. Default: 10"
		echo " "
		echo "	--cname:"
		echo "		Blosc compressor. Default: blosclz"
		echo " "
		echo "	--clevel:"
		echo "		Compression level (1 - 9). Default: 6"
		echo " "
		echo "	--shuffle:"
		echo "		Shuffle filter (NO_SHUFFLE, SHUFFLE, BIT_SHUFFLE). Default: SHUFFLE"
		echo " "
		echo "	--nthreads:"
		echo "		Set number of threads per task used."
		exit 0
		shift
		;;
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
		--cname)
		cname="$2"
		shift
		;;
		--clevel)
		clevel="$2"
		shift
		;;
		--shuffle)
		if [ "$2" == "NO_SHUFFLE" ]
		then
			shuffle=0
		elif [ "$2" == "SHUFFLE" ]
		then
			shuffle=1
		else
			shuffle=2
		fi
		shift
		;;
		--nthreads)
		nthreads="$2"
		shift
		;;
		*)
		echo USAGE: sbatch blosc_benchmark.sh
		echo " "
		echo "	-h | --help:"
		echo "		Shows this dialog."
		echo " "
		echo "	--in_name:"
		echo "		Input path to data."
		echo " "
		echo "	--format:"
		echo "		Data format (HACC, MPAS, EXALT)."
		echo " "
		echo "	--reps:"
		echo "		Number of repetitions for the benchmark."
		echo " "
		echo "	--cname:"
		echo "		Blosc compressor."
		echo " "
		echo "	--clevel:"
		echo "		Compression level (1 - 9)."
		echo " "
		echo "	--shuffle:"
		echo "		Shuffle filter (NO_SHUFFLE, SHUFFLE, BIT_SHUFFLE)."
		echo " "
		echo "	--nthreads:"
		echo "		Set number of threads per task used."
		exit 1
		shift
		;;
	esac
	shift
done

# create output folders if they do not exist yet
mkdir -p output
mkdir -p output/blosc_benchmark

# init result file
uuid=$(uuidgen)
uuid=${uuid,,}
result_name="output/blosc_benchmark/result_$uuid.csv"
if [ ! -f $result_name ]
then
	if [ "$format" == "HACC" ]
	then
		echo "in_name;format;reps;nranks;cname;clevel;shuffle;nthreads;osize;csize;ratio;x_ratio;y_ratio;z_ratio;vx_ratio;vy_ratio;vz_ratio;max_c_init;min_c_init;avg_c_init;max_c_exe;min_c_exe;avg_c_exe;max_c_clean;min_c_clean;avg_c_clean;max_d_init;min_d_init;avg_d_init;max_d_exe;min_d_exe;avg_d_exe;max_d_clean;min_d_clean;avg_d_clean" > $result_name
	fi
fi

# logging information
date
hostname
echo ---- BLOSC BENCHMARK ----
echo in_name:		$in_name
echo result_name:	$result_name
echo format:		$format
echo repetitions:	$reps
echo cname:		$cname
echo clevel:		$clevel
echo shuffle:		$shuffle
echo nthreads:		$nthreads
echo -------------------------
echo " "

# load modules
source ../env/bash.darwin.gcc-openmp

# load blosc library
export LD_LIBRARY_PATH=$(eval pwd)/../lib/3rdparty/blosc:$LD_LIBRARY_PATH

# openmp settings
export OMP_NUM_THREADS=$nthreads

# call program with the given parameters
mpirun ../build/bin/blosc_benchmark $in_name $result_name $format $reps $nthreads $clevel $shuffle $cname

#!/bin/bash
#SBATCH --output logs/gio_benchmark.%j.out
#SBATCH --error logs/gio_benchmark.%j.err
#SBATCH --partition=scaling

# init parameters
in_name=""
result_name=""
format=""
reps=10
num_buckets=1024
blocksize=2048
nthreads=1

# evaluate command line parameters
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		-h|--help)
		echo "USAGE: sbatch gio_benchmark.sh"
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
		echo "	--num_buckets:"
		echo "		Number of buckets used for the dynamic histogram. Default: 1024"
		echo " "
		echo "	--blocksize:"
		echo "		Block size of the data chunks to be compressed. Default: 2048"
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
		--num_buckets)
		num_buckets="$2"
		shift
		;;
		--blocksize)
		blocksize="$2"
		shift
		;;
		--nthreads)
		nthreads="$2"
		shift
		;;
		*)
		echo "USAGE: sbatch gio_benchmark.sh"
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
		echo "	--num_buckets:"
		echo "		Number of buckets used for the dynamic histogram. Default: 1024"
		echo " "
		echo "	--blocksize:"
		echo "		Block size of the data chunks to be compressed. Default: 2048"
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
mkdir -p output/gio_benchmark

# init result file
uuid=$(uuidgen)
uuid=${uuid,,}
result_name="output/gio_benchmark/result_$uuid.csv"
if [ ! -f $result_name ]
then
	if [ "$format" == "HACC" ]
	then
		echo "in_name;format;repetitions;number of ranks;num_buckets;blocksize;original size [B];compressed size [B];total compression ratio;x compression ratio;y compression ratio;z compression ratio;vx compression ratio;vy compression ratio;vz compression ratio;max compression init speed [GB/s];min compression init speed [GB/s];avg compression init speed [GB/s];max compression execute speed[GB/s];min compression execute speed [GB/s];avg compression execute speed [GB/s];max compression cleanup speed [GB/s];min compression cleanup speed [GB/s];avg compression cleanup speed [GB/s];max decompression init speed [GB/s];min decompression init speed [GB/s];avg decompression init speed [GB/s];max decompression execute speed [GB/s];min decompression execute speed [GB/s];avg decompression execute speed [GB/s];max decompression cleanup speed [GB/s];min decompression cleanup speed [GB/s];avg decompression cleanup speed [GB/s]" > $result_name
	fi
fi

# logging information
date
hostname
echo ---- GIO BENCHMARK ----
echo in_name:		$in_name
echo result_name:	$result_name
echo format:		$format
echo repetitions:	$reps
echo num_buckets:	$num_buckets
echo blocksize:		$blocksize
echo nthreads:		$nthreads
echo -----------------------
echo " "

# load modules
source ../env/bash.darwin.gcc-openmp

# openmp settings
export OMP_NUM_THREADS=$nthreads

# call program with the given parameters
mpirun ../bin/gio_benchmark $in_name $result_name $format $reps $num_buckets $blocksize

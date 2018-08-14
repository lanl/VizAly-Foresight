#!/bin/bash
#SBATCH --output logs/lzma_benchmark.%j.%t.out
#SBATCH --error logs/lzma_benchmark.%j.%t.err
#SBATCH --partition=scaling

# init parameters
in_name=""
result_name=""
format=""
reps=10
preset="9e"
nthreads=1

# evaluate command line parameters
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		-h|--help)
		echo "USAGE: sbatch lzma_benchmark.sh"
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
		echo "	--preset:"
		echo "		Compression level preset (1 - 9)[e]. Default: 9e"
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
		--preset)
		preset="$2"
		shift
		;;
		--nthreads)
		nthreads="$2"
		shift
		;;
		*)
		echo "USAGE: sbatch lzma_benchmark.sh"
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
		echo "	--preset:"
		echo "		Compression level preset (1 - 9)[e]. Default: 9e"
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
mkdir -p output/lzma_benchmark

# init result file
uuid=$(uuidgen)
uuid=${uuid,,}
result_name="output/lzma_benchmark/result_$uuid.csv"
if [ ! -f $result_name ]
then
	if [ "$format" == "HACC" ]
	then
		echo "in_name;format;reps;nranks;preset;nthreads;osize;csize;ratio;x_ratio;y_ratio;z_ratio;vx_ratio;vy_ratio;vz_ratio;max_c_init;min_c_init;avg_c_init;max_c_exe;min_c_exe;avg_c_exe;max_c_clean;min_c_clean;avg_c_clean;max_d_init;min_d_init;avg_d_init;max_d_exe;min_d_exe;avg_d_exe;max_d_clean;min_d_clean;avg_d_clean" > $result_name
	fi
fi

# logging information
date
hostname
echo ---- LZMA BENCHMARK ----
echo in_name:		$in_name
echo result_name:	$result_name
echo format:		$format
echo repetitions:	$reps
echo preset:		$preset
echo nthreads:		$nthreads
echo ------------------------
echo " "

# load modules
source ../env/bash.darwin.gcc-openmp

# load lzma library
export LD_LIBRARY_PATH=$(eval pwd)/../lib/3rdparty/lzma:$LD_LIBRARY_PATH

# openmp settings
export OMP_NUM_THREADS=$nthreads

# call program with the given parameters
mpirun ../bin/lzma_benchmark $in_name $result_name $format $reps $preset $nthreads

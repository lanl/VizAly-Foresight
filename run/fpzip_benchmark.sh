#!/bin/bash
#SBATCH --output logs/fpzip_benchmark.%j.%t.out
#SBATCH --error logs/fpzip_benchmark.%j.%t.err
#SBATCH --partition=scaling

# init parameters
in_name=""
out_name=""
result_name=""
format=""
reps=10
precision=32
nthreads=1

# evaluate command line parameters
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		-h|--help)
		echo "USAGE: sbatch fpzip_benchmark.sh"
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
		echo "	--precision:"
		echo "		Number of bits to remain for each value."
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
		--precision)
		precision="$2"
		shift
		;;
		--nthreads)
		nthreads="$2"
		shift
		;;
		*)
		echo "USAGE: sbatch fpzip_benchmark.sh"
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
		echo "	--precision:"
		echo "		Number of bits to remain for each value."
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
mkdir -p output/fpzip_benchmark

# init result file
uuid=$(uuidgen)
uuid=${uuid,,}
result_name="output/fpzip_benchmark/result_$uuid.csv"
if [ ! -f $result_name ]
then
	if [ "$format" == "HACC" ]
	then
		echo "in_name;format;reps;nranks;precision;nthreads;osize;csize;ratio;x_ratio;y_ratio;z_ratio;vx_ratio;vy_ratio;vz_ratio;x_abs_err;y_abs_err;z_abs_err;vx_abs_err;vy_abs_err;vz_abs_err;x_rel_err;y_rel_err;z_rel_err;vx_rel_err;vy_rel_err;vz_rel_err;x_snr;y_snr;z_snr;vx_snr;vy_snr;vz_snr;max_c_init;min_c_init;avg_c_init;max_c_exe;min_c_exe;avg_c_exe;max_c_clean;min_c_clean;avg_c_clean;max_d_init;min_d_init;avg_d_init;max_d_exe;min_d_exe;avg_d_exe;max_d_clean;min_d_clean;avg_d_clean" > $result_name
	fi
fi

# generate out_name
mkdir -p output/fpzip_benchmark/data
out_name="output/fpzip_benchmark/data/output_$precision.gio"

# logging information
date
hostname
echo ---- FPZIP BENCHMARK ----
echo in_name:		$in_name
echo out_name:		$out_name
echo result_name:	$result_name
echo format:		$format
echo repetitions:	$reps
echo precision:		$precision
echo nthreads:		$nthreads
echo -------------------------
echo " "

# load modules
source ../env/bash.darwin.gcc-openmp

# openmp settings
export OMP_NUM_THREADS=$nthreads

# call program with the given parameters
mpirun ../bin/fpzip_benchmark $in_name $out_name $result_name $format $reps $precision $nthreads

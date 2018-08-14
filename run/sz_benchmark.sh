#!/bin/bash
#SBATCH --output logs/sz_benchmark.%j.out
#SBATCH --error logs/sz_benchmark.%j.err
#SBATCH --partition=scaling

# init parameters
in_name=""
out_name=""
result_name=""
format="HACC"
reps=10
max_quant_intervals=2097152
quantization_intervals=0
data_endian_type=0
sol_id=101
layers=1
sample_distance=100
pred_threshold="0.99"
offset=0
sz_mode=1
gzip_mode=9
error_bound_mode=4
abs_err_bound=$(echo "10^-4" | bc -l)
rel_bound_ratio=$(echo "10-4" | bc -l)
pw_rel_bound_ratio=$(echo "10^-5" | bc -l)
segment_size=32
nthreads=1

# evaluate command line parameters
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		-h|--help)
		echo "USAGE: sbatch sz_benchmark.sh"
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
		echo "	--max_quant_intervals:"
		echo "		Maximum number of quantization intervals. Only valid if quantization_intervals=0. Default: 2097152"
		echo " "
		echo "	--quantization_intervals:"
		echo "		Number of quantization intervals has to be an even number. A value of 0 searches automatically for an optimized setting. Default: 0"
		echo " "
		echo "	--data_endian_type:"
		echo "		Endianess of the data (LITTLE_ENDIAN, BIG_ENDIAN). Default: LITTLE_ENDIAN"
		echo " "
		echo "	--layers:"
		echo "		Number of layers used in data prediction (1 - 3). Default: 1"
		echo " "
		echo "	--sample_distance:"
		echo "		Number of samples used to optimize the number of quantization intervals. Default: 100"
		echo " "
		echo "	--pred_threshold:"
		echo "		Threshold to determine the ratio of predictable data over all data. Default: 0.99"
		echo " "
		echo "	--offset:"
		echo "		Used to tune the compression ratio for hard-to-compress data (0 - 10). Default: 0"
		echo " "
		echo "	--sz_mode:"
		echo "		SZ mode (SZ_BEST_SPEED, SZ_DEFAULT_COMPRESSION, SZ_BEST_COMPRESSION). Default: SZ_BEST_COMPRESSION"
		echo " "
		echo "	--gzip_mode:"
		echo "		GZIP mode, only valid for SZ_BEST_COMPRESSION (Gzip_NO_COMPRESSION, Gzip_BEST_SPEED, Gzip_BEST_COMPRESSION, Gzip_DEFAULT_COMPRESSION). Default: Gzip_BEST_COMPRESSION"
		echo " "
		echo "	--error_bound_mode:"
		echo "		Error bound mode (ABS_AND_REL, ABS_OR_REL, ABS, REL, PW_REL). Default: PW_REL"
		echo " "
		echo "	--abs_err_bound:"
		echo "		Absolute error bound. Default: 10^-4"
		echo " "
		echo "	--rel_bound_ratio:"
		echo "		Relative bound ratio. Default: 10^-4"
		echo " "
		echo "	--pw_rel_bound_ratio:"
		echo "		Point-wise relative bound ratio. Default: 10^-4"
		echo " "
		echo "	--segment_size:"
		echo "		Point-wise relative bound ratio segment size. Default: 32"
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
		--max_quant_intervals)
		max_quant_intervals="$2"
		shift
		;;
		--quantization_intervals)
		quantization_intervals="$2"
		shift
		;;
		--data_endian_type)
		if [ "$2" == "LITTLE_ENDIAN" ]
		then
			data_endian_type=0
		else
			data_endian_type=1
		fi
		shift
		;;
		--layers)
		layers="$2"
		shift
		;;
		--sample_distance)
		sample_distance="$2"
		shift
		;;
		--pred_threshold)
		pred_threshold="$2"
		shift
		;;
		--offset)
		offset="$2"
		shift
		;;
		--sz_mode)
		if [ "$2" == "SZ_BEST_SPEED" ]
		then
			sz_mode=0
		elif [ "$2" == "SZ_DEFAULT_COMPRESSION" ]
		then
			sz_mode=2
		else
			sz_mode=1
		fi
		shift
		;;
		--gzip_mode)
		if [ "$2" == "Gzip_NO_COMPRESSION" ]
		then
			gzip_mode=0
		elif [ "$2" == "Gzip_BEST_SPEED" ]
		then
			gzip_mode=1
		elif [ "$2" == "Gzip_DEFAULT_COMPRESSION" ]
		then
			gzip_mode=-1
		else
			gzip_mode=9
		fi
		shift
		;;
		--error_bound_mode)
		if [ "$2" == "ABS_AND_REL" ]
		then
			error_bound_mode=2
		elif [ "$2" == "ABS_OR_REL" ]
		then
			error_bound_mode=3
		elif [ "$2" == "REL" ]
		then
			error_bound_mode=1
		elif [ "$2" == "PW_REL" ]
		then
			error_bound_mode=4
		else
			error_bound_mode=0
		fi
		shift
		;;
		--abs_err_bound)
		abs_err_bound=$(echo "$2" | bc -l)
		shift
		;;
		--rel_bound_ratio)
		rel_bound_ratio=$(echo "$2" | bc -l)
		shift
		;;
		--pw_rel_bound_ratio)
		pw_rel_bound_ratio=$(echo "$2" | bc -l)
		shift
		;;
		--segment_size)
		segment_size="$2"
		shift
		;;
		--nthreads)
		nthreads="$2"
		shift
		;;
		*)
		echo "$2 option does not exist!"
		echo " "
		echo " "
		echo "USAGE: sbatch sz_benchmark.sh"
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
		echo "	--max_quant_intervals:"
		echo "		Maximum number of quantization intervals. Only valid if quantization_intervals=0. Default: 2097152"
		echo " "
		echo "	--quantization_intervals:"
		echo "		Number of quantization intervals has to be an even number. A value of 0 searches automatically for an optimized setting. Default: 0"
		echo " "
		echo "	--data_endian_type:"
		echo "		Endianess of the data (LITTLE_ENDIAN, BIG_ENDIAN). Default: LITTLE_ENDIAN"
		echo " "
		echo "	--layers:"
		echo "		Number of layers used in data prediction (1 - 3). Default: 1"
		echo " "
		echo "	--sample_distance:"
		echo "		Number of samples used to optimize the number of quantization intervals. Default: 100"
		echo " "
		echo "	--pred_threshold:"
		echo "		Threshold to determine the ratio of predictable data over all data. Default: 0.99"
		echo " "
		echo "	--offset:"
		echo "		Used to tune the compression ratio for hard-to-compress data (0 - 10). Default: 0"
		echo " "
		echo "	--sz_mode:"
		echo "		SZ mode (SZ_BEST_SPEED, SZ_DEFAULT_COMPRESSION, SZ_BEST_COMPRESSION). Default: SZ_BEST_COMPRESSION"
		echo " "
		echo "	--gzip_mode:"
		echo "		GZIP mode, only valid for SZ_BEST_COMPRESSION (Gzip_NO_COMPRESSION, Gzip_BEST_SPEED, Gzip_BEST_COMPRESSION, Gzip_DEFAULT_COMPRESSION). Default: Gzip_BEST_COMPRESSION"
		echo " "
		echo "	--error_bound_mode:"
		echo "		Error bound mode (ABS_AND_REL, ABS_OR_REL, ABS, REL, PW_REL). Default: PW_REL"
		echo " "
		echo "	--abs_err_bound:"
		echo "		Absolute error bound. Default: 10^-4"
		echo " "
		echo "	--rel_bound_ratio:"
		echo "		Relative bound ratio. Default: 10^-4"
		echo " "
		echo "	--pw_rel_bound_ratio:"
		echo "		Point-wise relative bound ratio. Default: 10^-4"
		echo " "
		echo "	--segment_size:"
		echo "		Point-wise relative bound ratio segment size. Default: 32"
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
mkdir -p output/sz_benchmark

# init result file
uuid=$(uuidgen)
uuid=${uuid,,}
result_name="output/sz_benchmark/result_$uuid.csv"
if [ ! -f $result_name ]
then
	if [ "$format" == "HACC" ]
	then
		echo "in_name;format;repetitions;nranks;max_quant_intervals;quantization_intervals;data_endian_type;sol_id;layers;sample_distance;pred_threshold;offset;sz_mode;gzip_mode;error_bound_mode;abs_err_bound;rel_bound_ratio;pw_rel_bound_ratio;segment_size;nthreads;osize;csize;ratio;x_ratio;y_ratio;z_ratio;vx_ratio;vy_ratio;vz_ratio;x_abs_err;y_abs_err;z_abs_err;vx_abs_err;vy_abs_err;vz_abs_err;x_rel_err;y_rel_err;z_rel_err;vx_rel_err;vy_rel_err;vz_rel_err;x_snr;y_snr;z_snr;vx_snr;vy_snr;vz_snr;max_c_init;min_c_init;avg_c_init;max_c_exe;min_c_exe;avg_c_exe;max_c_clean;min_c_clean;avg_c_clean;max_d_init;min_d_init;avg_d_init;max_d_exe;min_d_exe;avg_d_exe;max_d_clean;min_d_clean;avg_d_clean" > $result_name
	fi
fi

# generate out_name
mkdir -p output/sz_benchmark/data
out_name="output/sz_benchmark/data/output_$uuid.csv"

# logging information
date
hostname
echo ---- SZ BENCHMARK ----
echo in_name:			$in_name
echo out_name:			$out_name
echo result_name:		$result_name
echo format:			$format
echo repetitions:		$reps
echo max_quant_intervals:	$max_quant_intervals
echo quantization_intervals:	$quantization_intervals
echo data_endian_type:		$data_endian_type
echo sol_id:			$sol_id
echo layers:			$layers
echo sample_distance:		$sample_distance
echo pred_threshold:		$pred_threshold
echo offset:			$offset
echo sz_mode:			$sz_mode
echo gzip_mode:			$gzip_mode
echo error_bound_mode:		$error_bound_mode
echo abs_err_bound:		$abs_err_bound
echo rel_bound_ratio:		$rel_bound_ratio
echo pw_rel_bound_ratio:	$pw_rel_bound_ratio
echo segment_size:		$segment_size
echo nthreads:			$nthreads
echo ---------------------
echo " "

# load modules
source ../env/bash.darwin.gcc-openmp

# load sz library
export LD_LIBRARY_PATH=$(eval pwd)/../lib/3rdparty/sz:$LD_LIBRARY_PATH

# openmp settings
export OMP_NUM_THREADS=$nthreads

# call program with the given parameters
mpirun ../bin/sz_benchmark $in_name $out_name $result_name $format $reps $max_quant_intervals $quantization_intervals $data_endian_type $sol_id $layers $sample_distance $pred_threshold $offset $sz_mode $gzip_mode $error_bound_mode $abs_err_bound $rel_bound_ratio $pw_rel_bound_ratio $segment_size $nthreads -n

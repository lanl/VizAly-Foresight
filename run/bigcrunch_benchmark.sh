#!/bin/bash
#SBATCH --output logs/bigcrunch_benchmark.%j.%t.out
#SBATCH --error logs/bigcrunch_benchmark.%j.%t.err
#SBATCH --partition=scaling

# init parameters
in_name=""
out_name=""
result_name=""
reps=10
nthreads=1
error=-3
tolerance=0

# evaluate command line parameters
while [[ $# -gt 0 ]]
do
	key="$1"
	case $key in
		-h|--help)
		echo "USAGE: sbatch bigcrunch_benchmark.sh"
		echo " "
		echo "	-h | --help:"
		echo "		Shows this dialog."
		echo " "
		echo "	--in_name:"
		echo "		Input path to data."
		echo " "
		echo "	--reps:"
		echo "		Number of repetitions for the benchmark. Default: 10"
		echo " "
		echo "	--error:"
		echo "		Target relative error. Default: -3"
		echo " "
		echo "	--tolerance:"
		echo "		Exponent at which relative error is switched with absolute error. Default: 0"
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
		--reps)
		reps="$2"
		shift
		;;
		--error)
		error="$2"
		shift
		;;
		--tolerance)
		tolerance="$2"
		shift
		;;
		--nthreads)
		nthreads="$2"
		shift
		;;
		*)
		echo "USAGE: sbatch bigcrunch_benchmark.sh"
		echo " "
		echo "	-h | --help:"
		echo "		Shows this dialog."
		echo " "
		echo "	--in_name:"
		echo "		Input path to data."
		echo " "
		echo "	--reps:"
		echo "		Number of repetitions for the benchmark. Default: 10"
		echo " "
		echo "	--error:"
		echo "		Target relative error. Default: -3"
		echo " "
		echo "	--tolerance:"
		echo "		Exponent at which relative error is switched with absolute error. Default: 0"
		echo " "
		echo "	--nthreads:"
		echo "		Set number of threads per task used."
		exit 1
		shift
		;;
	esac
	shift
done

# create output folder if they do not exist yet
mkdir -p output
mkdir -p output/bigcrunch_benchmark

# init result file
uuid=$(uuidgen)
uuid=${uuid,,}
result_name="output/bigcrunch_benchmark/result_$uuid.csv"
if [ ! -f $result_name ]
then
	echo "in_name,reps,nranks,error,tolerance,nthreads,osize,csize,ratio,x_ratio,y_ratio,z_ratio,vx_ratio,vy_ratio,vz_ratio,x_rel_err,y_rel_err,z_rel_err,vx_rel_err,vy_rel_err,vz_rel_err,max_c_init,min_c_init,avg_c_init,max_c_exe,min_c_exe,avg_c_exe,max_c_clean,min_c_clean,avg_c_clean,max_d_init,min_d_init,avg_d_init,max_d_exe,min_d_exe,avg_d_exe,max_d_clean,min_d_clean,avg_d_clean" > $result_name
fi

# generate out_name
#mkdir -p output/bigcrunch_benchmark/data
#out_name="output/bigcrunch_benchmark/data/output_"$error"_"$tolerance".gio"

# logging information
date
hostname
echo ---- BIGCRUNCH BENCHMARK ----
echo in_name:		$in_name
echo out_name:		$out_name
echo result_name:	$result_name
echo repetitions:	$reps
echo error:		$error
echo tolerance:		$tolerance
echo nthreads:		$nthreads
echo -----------------------------
echo " "

# load modules
source ../env/bash.darwin.gcc-openmp

# openmp settings
export OMP_NUM_THREADS=$nthreads

# call program with the given parameters
mpirun -n 1 ../build/bin/bigcrunch_benchmark $in_name $result_name $reps $nthreads $error $tolerance

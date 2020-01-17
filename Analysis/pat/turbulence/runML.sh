#!/bin/bash
 
#SBATCH -N 1
#SBATCH -p volta-x86
#SBATCH --job-name SZ_abs_0.1

pwd; hostname; date


# project variables
env="/usr/projects/ml_compression/VizAly-Foresight-CBench/scripts/VizAly-CBench.bash.darwin-autoencoder"

project_path="/usr/projects/ml_compression/VizAly-Foresight-CBench/Analysis/pat/turbulence/"
run_path="/usr/projects/ml_compression/runs/"

original_data="/usr/projects/ml_compression/data/scalarHIT_fields100.h5"
test_data="/usr/projects/ml_compression/SZ_abs_0.01___turbulence.h5"

output_dir="diagnosticsCAE_SZ_0.1_"
plot_results_dir="autoencoder_SZ_0.1_"
output_prefix="SZ_0.1__"


# Initialize environment
source ${env}
cd ${run_path}
mkdir ${output_dir}

# Run
python ${project_path}"turbulence_analytics/debugInference.py" ${original_data} ${test_data} ${output_dir}
python ${project_path}"turbulence_analytics/plot_all.py" ${output_dir} ${output_prefix} ${plot_results_dir}
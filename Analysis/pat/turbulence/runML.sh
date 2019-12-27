#!/bin/bash
 
#SBATCH -N 1
#SBATCH -p volta-x86

pwd; hostname; date


# project variables
env="/projects/ml_compression/VizAly-Foresight/scripts/VizAly-CBench.bash.darwin-autoencoder"

project_path="/projects/ml_compression/VizAly-Foresight/Analysis/pat/turbulence/"
run_path="/projects/ml_compression/VizAly-Foresight/Analysis/pat/turbulence/"

original_data="/projects/ml_compression/data/scalarHIT_fields100.h5"
test_data="/projects/ml_compression/data/raw50.h5"

output_dir="diagnosticsCAE1"
plot_results_dir="autoencoder_Oct_23"
output_prefix="test_"


# Initialize environment
source ${env}
cd ${run_path}
mkdir ${output_dir}

# Run
python ${project_path}"turbulence_analytics/debugInference.py" ${original_data} ${test_data $output_dir}
python ${project_path}"turbulence_analytics/plot_all.py" ${output_dir} ${output_prefix} ${plot_results_dir}
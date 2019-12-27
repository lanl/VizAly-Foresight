#!/bin/bash
 
#SBATCH --partition=scaling
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8


pwd; hostname; date


# project variables
project_path="/usr/projects/ml_compression/VizAly-Foresight-CBench/"
run_path="/usr/projects/ml_compression/runs/"

exec="build/CBench"

env="scripts/VizAly-CBench.bash.darwin"
output_dir="test-ml-turbulence"


# Initialize environment
source ${project_path}${env}
cd ${run_path}
mkdir ${output_dir}

# Run
json_path=${project_path}"inputs/raw/ml_turbulence.json"
echo ${json_path}
echo ${project_path}${exec}
mpirun ${project_path}${exec} ${json_path}
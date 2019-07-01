#!/bin/bash
#SBATCH -N 2
#SBATCH --ntasks-per-node 4
#SBATCH -p general

projectPath=$(pwd)

# load modules
source $projectPath/scripts/VizAly-CBench.bash.darwin

# go to folder 
cd $projectPath/build

# Run:
mpirun $projectPath/build/CBench $projectPath/inputs/hacc/hacc_cbench_test.json
#!/bin/bash

#SBATCH --nodes 8
#SBATCH --ntasks-per-node 2
#SBATCH --partition scaling
#SBATCH --job-name merge

BUILD_DIR="/projects/exasky/hoby-projects/cbench/build"

# load modules
source "/home/hoby/.bashrc"

# go to run folder
cd "${BUILD_DIR}"

# run
mpirun -np 16 ./merger ../inputs/hacc/hacc_analysis_merge.json

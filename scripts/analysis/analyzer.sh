#!/bin/bash

#SBATCH --nodes 8
#SBATCH --ntasks-per-node 16
#SBATCH --partition scaling
#SBATCH --job-name extract

BUILD_DIR="/projects/exasky/hoby-projects/cbench/build"

# load modules
source "/home/hoby/.bashrc"

# go to run folder
cd "${BUILD_DIR}"

# run
mpirun -np 128 ./analyzer ../inputs/hacc/hacc_analysis_non-halo.json

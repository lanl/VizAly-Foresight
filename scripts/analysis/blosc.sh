#!/bin/bash

#SBATCH --nodes 4
#SBATCH --ntasks-per-node 4
#SBATCH --partition scaling
#SBATCH --job-name halos

HACC="/projects/exasky/HACC"
BUILD_DIR="/projects/exasky/hoby-projects/cbench/build"

# load modules
source "${HACC}.darwin_setup"

# go to folder 
cd ${BUILD_DIR}

# run
mpirun -np 16 ./CBench ../inputs/hacc/hacc_analysis_blosc_halo.json

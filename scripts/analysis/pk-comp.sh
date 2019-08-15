#!/bin/bash

#SBATCH --nodes 4
#SBATCH --ntasks-per-node 2
#SBATCH --partition scaling
#SBATCH --job-name spectrum

HACC="/projects/exasky/HACC"
POWER_SPECTRUM="${HACC}/trunk/Darwin/mpi/bin/hacc_pk_gio_auto"
PARTICLES_DATA="/projects/exasky/hoby-projects/cbench/build/data/all_combined_m000-499.particles"
OUTPUT_DATA="/projects/exasky/hoby-projects/cbench/build/analysis-decompressed-pk-499.dat"
TIMESTEP=499

# load modules
source "${HACC}.darwin_setup"

# go to run folder
cd "${HACC}/run"

# run
mpirun ${POWER_SPECTRUM} inputs/indat.params -n ${PARTICLES_DATA} ${OUTPUT_DATA} ${TIMESTEP}

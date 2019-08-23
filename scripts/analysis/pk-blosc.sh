#!/bin/bash

#SBATCH --nodes 4
#SBATCH --ntasks-per-node 4
#SBATCH --partition scaling
#SBATCH --job-name spectrum

HACC="/projects/exasky/HACC"
POWER_SPECTRUM="${HACC}/trunk/Darwin/mpi/bin/hacc_pk_gio_auto"
PARTICLES_DATA="/projects/exasky/hoby-projects/cbench/build/data/blosc_combined_m000-499.particles"
OUTPUT_DATA="/projects/exasky/hoby-projects/cbench/build/data/analysis-spectrum-full-blosc.dat"
TIMESTEP=499

# load modules
source "${HACC}.darwin_setup"

# go to run folder
cd "${HACC}/run"

# run
mpirun ${POWER_SPECTRUM} inputs/indat.params -n ${PARTICLES_DATA} ${OUTPUT_DATA} ${TIMESTEP}

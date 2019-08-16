#!/bin/bash

#SBATCH --nodes 4
#SBATCH --ntasks-per-node 4
#SBATCH --partition scaling
#SBATCH --job-name spectrum

HACC="/projects/exasky/HACC"
POWER_SPECTRUM="${HACC}/trunk/Darwin/mpi/bin/hacc_pk_gio_auto"
#PARTICLES_DATA="/projects/exasky/data/Argonne_L360_HACC001/STEP499/m000.full.mpicosmo.499"
#PARTICLES_DATA="/projects/exasky/hoby-projects/cbench/build/data/all_decompressed_sz__m000.full.mpicosmo.499"
PARTICLES_DATA="/projects/exasky/hoby-projects/cbench/build/data/all_decompressed_fpzip__m000.full.mpicosmo.499"
OUTPUT_DATA="/projects/exasky/hoby-projects/cbench/build/data/analysis-spectrum-full-fpzip.dat"
#OUTPUT_DATA="${HACC}/run/output/analysis/Pow/uncompressed.pk.499"
TIMESTEP=499

# load modules
source "${HACC}.darwin_setup"

# go to run folder
cd "${HACC}/run"

# run
mpirun ${POWER_SPECTRUM} inputs/indat.params -n ${PARTICLES_DATA} ${OUTPUT_DATA} ${TIMESTEP}

#!/bin/bash

#SBATCH --nodes 8
#SBATCH --ntasks-per-node 1
#SBATCH --partition skylake-gold
#SBATCH --job-name analysis

HACC="/projects/exasky/HACC"
CBENCH="/projects/exasky/hoby-projects/cbench"
BUILD="${CBENCH}/build"
POWER_SPECTRUM="${HACC}/trunk/Darwin/mpi/bin/hacc_pk_gio_auto"
PARTICLES_DATA="${BUILD}/data/data-combined-zip"
OUTPUT_DATA="${BUILD}/data/pk-combined-zip.dat"
TIMESTEP=499

# load modules
source "${HACC}.darwin_setup" && cd ${BUILD} &&

# compress non-halo particles dataset
mpirun -np 8 ./CBench ../inputs/hacc/hacc_analysis_compress_non-halo.json &&

# merge it with halo ones
mpirun -np 8 ./merger ../inputs/hacc/hacc_analysis_merge.json &&

# compute power spectrum eventually
cd "${HACC}/run" &&
mpirun ${POWER_SPECTRUM} inputs/indat.params -n ${PARTICLES_DATA} ${OUTPUT_DATA} ${TIMESTEP}
#!/bin/bash

#SBATCH --nodes 8
#SBATCH --ntasks-per-node 1
#SBATCH --partition skylake-gold
#SBATCH --job-name analysis

# enable or disable steps
EXTRACT_NON_HALOS=true
COMPRESS_NON_HALOS=false
MERGE_DATASETS=false
COMPUTE_POWER_SPECTRUM=false

# parameters
NRANKS=8
TIMESTEP=499
SUFFIX="sampled"
HACC="/projects/exasky/HACC"
BUILD="/projects/exasky/hoby-projects/cbench/build"
INPUT_JSON="../inputs/hacc/analysis_pipeline_${SUFFIX}.json"
POWER_SPECTRUM="${HACC}/trunk/Darwin/mpi/bin/hacc_pk_gio_auto"
PARTICLES_DATA="/projects/exasky_data/hoby/analysis/data-combined-zip-${SUFFIX}-0.90"
OUTPUT_DATA="/projects/exasky_data/hoby/analysis/pk-combined-zip-${SUFFIX}-0.90.dat"
STATUS=0

# extract non-halos and compute entropy if required
if ${EXTRACT_NON_HALOS}; then
  source "/home/hoby/.bashrc" && cd ${BUILD} &&
  mpirun -np ${NRANKS} ./analyzer ${INPUT_JSON} &&
  STATUS=$?
fi

# compress non-halo particles dataset
if ${COMPRESS_NON_HALOS} && [ ${STATUS} -eq 0 ]; then
  source "/home/hoby/.bashrc" && cd ${BUILD} &&
  mpirun -np ${NRANKS} ./CBench ${INPUT_JSON} &&
  STATUS=$?
fi

# merge it with halo ones
if ${MERGE_DATASETS} && [ ${STATUS} -eq 0 ]; then
  source "/home/hoby/.bashrc" && cd ${BUILD} &&
  mpirun -np ${NRANKS} ./merger ${INPUT_JSON} &&
  STATUS=$?
fi

# compute power spectrum eventually
if ${COMPUTE_POWER_SPECTRUM} && [ ${STATUS} -eq 0 ]; then
  source "${HACC}.darwin_setup" && cd "${HACC}/run" &&
  mpirun ${POWER_SPECTRUM} inputs/indat.params -n ${PARTICLES_DATA} ${OUTPUT_DATA} ${TIMESTEP}
fi
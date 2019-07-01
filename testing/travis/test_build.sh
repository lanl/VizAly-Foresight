#! /bin/bash

set -e

# get GCC compiler
source /etc/bashrc
source /src/env.sh

export PLATFORM="travis"

# load OpenMPI
module load mpi

# build
cd /src/VizAly-CBench
projectPath=$(pwd)

source buildDependencies.sh # build dependencies
source build.sh	-travis    	# build the code

# run example
mpirun -np 4 --allow-run-as-root ./CBench ../testing/scripts/hacc_input.json
mpirun -np 4 --allow-run-as-root ./CBench ../testing/scripts/nyx_input.json

# view output
cat metrics_HACC_Travis_.csv
cat metrics_NYX_Travis_.csv


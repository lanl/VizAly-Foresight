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
source build.sh    			# build the code

# run example
mpirun -np 4 --allow-run-as-root ./CBench ../testing/scripts/hacc_input.json
mpirun -np 4 --allow-run-as-root ./CBench ../testing/scripts/nyx_input.json

# view output
cat metrics_HACC_Travis_.csv
echo ""
cat metrics_NYX_Travis_.csv

# install Python dependencies
cd ..
pip3 install -r testing/travis/requirements.txt


# test HACC workflow generator and executables
cd /src/VizAly-CBench/Analysis
python3 -m pat.nyx.workflow --input-file ../testing/scripts/NYX_workflow.json --preview

#python3 -m pat.hacc.pat_hacc --input-file ../inputs/hacc/hacc_wflow.json --preview
#python3 -m pat.hacc.pat_hacc_query.py --help
#python3 -m pat.hacc.pat_cinema.py --help

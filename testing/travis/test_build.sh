#! /bin/bash

set -e

# get GCC compiler
source /src/env.sh

# build
cd /src/VizAly-CBench
projectPath=$(pwd)

source buildDependencies.sh	# build dependencies
source build.sh				# build the code

# run example
mpirun -np 4 --allow-run-as-root ./CBench ../testing/hacc_input.json
#mpirun -np 4 --allow-run-as-root ./CBench ../testing/nyx_input.json

# Compare output with previous run
python ../testing/travis/test_output.py

# view output
cat metrics_HACC_all_

# 
